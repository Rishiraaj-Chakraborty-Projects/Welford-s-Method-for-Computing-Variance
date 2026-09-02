#include "welford/accumulator.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

static int check(const char* name, bool ok) {
    std::cout << "  [ " << (ok ? "PASS" : "FAIL") << " ] " << name << std::endl;
    return ok ? 0 : 1;
}

int main() {
    int failures = 0;
    std::cout << "Running Welford Accumulator Unit Tests...\n" << std::endl;

    // Test 1: Textbook Sample
    {
        std::vector<double> sample = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
        welford::WelfordAccumulator acc;
        for (double v : sample) acc.update(v);
        
        failures += check("Textbook mean calculation", std::abs(acc.mean() - 5.0) < 1e-12);
        failures += check("Textbook sample variance calculation", std::abs(acc.variance() - (32.0 / 7.0)) < 1e-12);
    }

    // Test 2: Catastrophic Cancellation Series
    {
        const auto series = welford::generate_cancellation_series(100000, 1e8);
        double ref_mean = 0.0, ref_var = 0.0;
        welford::compute_reference_stats(series, ref_mean, ref_var);

        welford::WelfordAccumulator welford_acc;
        welford::NaiveAccumulator naive_acc;
        for (double v : series) {
            welford_acc.update(v);
            naive_acc.update(v);
        }

        const double welford_rel_err = std::abs(welford_acc.variance() - ref_var) / ref_var;
        const double naive_rel_err = std::abs(naive_acc.variance() - ref_var) / ref_var;

        failures += check("Welford within 1e-6 relative error on high-mean dataset", welford_rel_err < 1e-6);
        failures += check("Naive algorithm demonstrates catastrophic cancellation error (> 0.5 relative error)", naive_rel_err > 0.5);
    }

    // Test 3: Higher Moments (Uniform Distribution Skew & Kurtosis)
    {
        welford::MomentsAccumulator moments_acc;
        for (int i = 0; i < 100000; ++i) {
            const double val = static_cast<double>(i % 1000) / 999.0;
            moments_acc.update(val);
        }
        failures += check("Uniform distribution skewness near zero", std::abs(moments_acc.skewness()) < 0.02);
        failures += check("Uniform distribution excess kurtosis near -1.2", std::abs(moments_acc.excess_kurtosis() + 1.2) < 0.05);
    }

    // Test 4: Bivariate Linear Pair
    {
        welford::BivariateAccumulator biv_acc;
        for (int i = 0; i < 1000; ++i) {
            const double v = static_cast<double>(i);
            biv_acc.update(v, 3.0 * v + 1.0);
        }
        failures += check("Bivariate linear pair yields beta = 3.0", std::abs(biv_acc.beta() - 3.0) < 1e-9);
        failures += check("Bivariate linear pair yields correlation = 1.0", std::abs(biv_acc.correlation() - 1.0) < 1e-9);
    }

    // Test 5: Shift Invariance
    {
        welford::WelfordAccumulator acc_base, acc_shifted;
        for (int i = 0; i < 1000; ++i) {
            acc_base.update(static_cast<double>(i));
            acc_shifted.update(static_cast<double>(i) + 1e6);
        }
        const double rel_diff = std::abs(acc_base.variance() - acc_shifted.variance()) / acc_base.variance();
        failures += check("Variance calculation is shift-invariant", rel_diff < 1e-9);
    }

    // Test 6: Single Observation Edge Case
    {
        welford::WelfordAccumulator acc;
        acc.update(42.0);
        failures += check("Single observation reports zero sample variance", acc.variance() == 0.0);
    }

    std::cout << "\nSummary: " << (failures == 0 ? "ALL TESTS PASSED" : "TEST FAILURES DETECTED") << std::endl;
    return failures == 0 ? 0 : 1;
}
