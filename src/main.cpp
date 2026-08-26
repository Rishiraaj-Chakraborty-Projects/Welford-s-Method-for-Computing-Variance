#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

struct Naive {
    double sum = 0, sumsq = 0; long n = 0;
    void push(double x) { sum += x; sumsq += x * x; ++n; }
    double mean() const { return n ? sum / n : 0.0; }
    double var() const { return sumsq / n - (sum / n) * (sum / n); }
};

struct Welford {
    double mean_ = 0, m2 = 0; long n = 0;
    void push(double x) {
        ++n;
        const double d = x - mean_;
        mean_ += d / n;
        m2 += d * (x - mean_);
    }
    double mean() const { return mean_; }
    double var() const { return n > 1 ? m2 / (n - 1) : 0.0; }
    double sd() const { return std::sqrt(var()); }
};

// same update carried out to the fourth moment. the extra terms are what you
// need for skew and kurtosis without a second pass over the data.
struct Moments {
    double mean_ = 0, m2 = 0, m3 = 0, m4 = 0; long n = 0;
    void push(double x) {
        const long n1 = n;
        ++n;
        const double d = x - mean_;
        const double dn = d / n;
        const double dn2 = dn * dn;
        const double t = d * dn * n1;
        mean_ += dn;
        m4 += t * dn2 * (n * n - 3.0 * n + 3.0) + 6.0 * dn2 * m2 - 4.0 * dn * m3;
        m3 += t * dn * (n - 2.0) - 3.0 * dn * m2;
        m2 += t;
    }
    double mean() const { return mean_; }
    double var() const { return n > 1 ? m2 / (n - 1) : 0.0; }
    double skew() const { return n > 2 ? std::sqrt(double(n)) * m3 / std::pow(m2, 1.5) : 0.0; }
    double kurt() const { return n > 3 ? double(n) * m4 / (m2 * m2) - 3.0 : 0.0; }
};

// two variable version. the covariance term is the numerator of a hedge ratio,
// which is the reason this matters on a desk rather than in a textbook.
struct CoWelford {
    double mx = 0, my = 0, mxx = 0, myy = 0, mxy = 0; long n = 0;
    void push(double x, double y) {
        ++n;
        const double dx = x - mx, dy = y - my;
        mx += dx / n;
        my += dy / n;
        mxx += dx * (x - mx);
        myy += dy * (y - my);
        mxy += dx * (y - my);
    }
    double cov() const { return n > 1 ? mxy / (n - 1) : 0.0; }
    double corr() const {
        const double d = std::sqrt(mxx * myy);
        return d > 0 ? mxy / d : 0.0;
    }
    double beta() const { return mxx > 0 ? mxy / mxx : 0.0; }
};

// exponentially weighted, the estimator most risk systems actually run
struct Ewma {
    double lambda, mean_ = 0, var_ = 0; bool init = false;
    explicit Ewma(double l) : lambda(l) {}
    void push(double x) {
        if (!init) { mean_ = x; var_ = 0; init = true; return; }
        const double d = x - mean_;
        mean_ += (1.0 - lambda) * d;
        var_ = lambda * var_ + (1.0 - lambda) * d * d;
    }
    double var() const { return var_; }
    double sd() const { return std::sqrt(var_); }
};

static void reference(const std::vector<double>& x, double& mean, double& var) {
    double m = 0;
    for (double v : x) m += v;
    m /= x.size();
    double t = 0;
    for (double v : x) t += (v - m) * (v - m);
    mean = m;
    var = t / (x.size() - 1);
}

static std::vector<double> cancellation_series(long n, double level) {
    std::vector<double> x;
    x.reserve(n);
    for (long i = 0; i < n; ++i)
        x.push_back(level + ((i % 7) - 3) * 1e-3);
    return x;
}

static int demo(long n, double level) {
    const std::vector<double> x = cancellation_series(n, level);

    Naive a; Welford b; Moments c;
    for (double v : x) { a.push(v); b.push(v); c.push(v); }

    double rm, rv;
    reference(x, rm, rv);

    std::printf("%ld observations of size %.0e with a spread of 1e-3 around it\n\n", n, level);
    std::printf("  %-24s %20s %14s\n", "method", "variance", "rel error");
    std::printf("  %-24s %20.10e %14s\n", "two pass (reference)", rv, "-");
    std::printf("  %-24s %20.10e %14.3e\n", "naive sum of squares", a.var(),
                std::abs(a.var() - rv) / rv);
    std::printf("  %-24s %20.10e %14.3e\n", "Welford", b.var(),
                std::abs(b.var() - rv) / rv);
    std::printf("  %-24s %20.10e %14.3e\n", "Welford with moments", c.var(),
                std::abs(c.var() - rv) / rv);

    std::printf("\n  the naive figure is wrong by a factor of %.0e. it is not noisy,\n",
                std::abs(a.var() - rv) / rv);
    std::printf("  it is meaningless. sumsq and sum*sum agree to about 16 digits\n");
    std::printf("  and the answer lives in the digits below that.\n");

    std::printf("\nhigher moments, same single pass\n\n");
    std::printf("  mean      %.10e\n", c.mean());
    std::printf("  variance  %.10e\n", c.var());
    std::printf("  skewness  %+.6f\n", c.skew());
    std::printf("  kurtosis  %+.6f (excess)\n", c.kurt());

    std::printf("\nonline covariance on a synthetic hedge pair\n\n");
    CoWelford cw;
    unsigned long s = 12345;
    auto u = [&s] {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return double((s >> 11) & ((1ULL << 53) - 1)) / double(1ULL << 53);
    };
    auto z = [&u] {
        const double a1 = u() + 1e-12, a2 = u();
        return std::sqrt(-2.0 * std::log(a1)) * std::cos(6.283185307179586 * a2);
    };
    const double true_beta = 1.35, resid_sd = 0.4;
    for (int i = 0; i < 200000; ++i) {
        const double f = z();
        cw.push(f, true_beta * f + resid_sd * z());
    }
    std::printf("  true beta       %.4f\n", true_beta);
    std::printf("  estimated beta  %.4f\n", cw.beta());
    std::printf("  correlation     %.4f\n", cw.corr());
    std::printf("  implied 1-rho^2 %.4f  (residual share of variance)\n",
                1.0 - cw.corr() * cw.corr());

    std::printf("\nexponentially weighted variance on a regime change\n\n");
    Ewma e94(0.94), e97(0.97);
    std::printf("  %-6s %14s %14s\n", "step", "lambda 0.94", "lambda 0.97");
    for (int i = 0; i < 400; ++i) {
        const double sd = i < 200 ? 0.01 : 0.04;
        const double r = sd * z();
        e94.push(r); e97.push(r);
        if (i == 199 || i == 205 || i == 220 || i == 260 || i == 399)
            std::printf("  %-6d %14.6f %14.6f\n", i + 1, e94.sd(), e97.sd());
    }
    std::printf("\n  true sd steps from 0.01 to 0.04 at step 200. the faster decay\n");
    std::printf("  tracks the break sooner and is noisier once it gets there.\n");
    return 0;
}

static int check(const char* name, bool ok) {
    std::printf("  %-46s %s\n", name, ok ? "ok" : "FAIL");
    return ok ? 0 : 1;
}

static int selftest() {
    int bad = 0;
    std::printf("self test\n\n");

    {
        std::vector<double> x = {2, 4, 4, 4, 5, 5, 7, 9};
        Welford w; for (double v : x) w.push(v);
        bad += check("mean of the textbook sample", std::abs(w.mean() - 5.0) < 1e-12);
        bad += check("sample variance of the same", std::abs(w.var() - 32.0 / 7.0) < 1e-12);
    }
    {
        const std::vector<double> x = cancellation_series(100000, 1e8);
        double rm, rv; reference(x, rm, rv);
        Welford w; Naive nv;
        for (double v : x) { w.push(v); nv.push(v); }
        bad += check("Welford within 1e-6 relative on hard data",
                     std::abs(w.var() - rv) / rv < 1e-6);
        bad += check("naive blows up by more than 1e6 relative",
                     std::abs(nv.var() - rv) / rv > 1e6);
    }
    {
        Moments m;
        for (int i = 0; i < 100000; ++i) {
            const double t = (i % 1000) / 999.0;
            m.push(t);
        }
        bad += check("uniform skewness near zero", std::abs(m.skew()) < 0.02);
        bad += check("uniform excess kurtosis near -1.2",
                     std::abs(m.kurt() + 1.2) < 0.05);
    }
    {
        CoWelford cw;
        for (int i = 0; i < 1000; ++i) { const double v = i; cw.push(v, 3.0 * v + 1.0); }
        bad += check("perfect linear pair gives beta 3", std::abs(cw.beta() - 3.0) < 1e-9);
        bad += check("perfect linear pair gives corr 1", std::abs(cw.corr() - 1.0) < 1e-9);
    }
    {
        Welford a, b;
        for (int i = 0; i < 1000; ++i) { a.push(i); b.push(i + 1e6); }
        bad += check("variance is shift invariant",
                     std::abs(a.var() - b.var()) / a.var() < 1e-9);
    }
    {
        Welford w; w.push(1.0);
        bad += check("single observation reports zero variance", w.var() == 0.0);
    }

    std::printf("\n%s\n", bad ? "FAILURES PRESENT" : "all checks passed");
    return bad ? 1 : 0;
}

int main(int argc, char** argv) {
    long n = 100000;
    double level = 1e8;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--test")) return selftest();
        if (!std::strcmp(argv[i], "--help")) {
            std::printf("usage: welford [--n N] [--level L] [--test]\n");
            return 0;
        }
        if (!std::strcmp(argv[i], "--n") && i + 1 < argc) n = std::atol(argv[++i]);
        else if (!std::strcmp(argv[i], "--level") && i + 1 < argc) level = std::atof(argv[++i]);
    }
    return demo(n, level);
}
