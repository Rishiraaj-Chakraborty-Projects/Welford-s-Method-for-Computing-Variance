#ifndef WELFORD_ACCUMULATOR_HPP
#define WELFORD_ACCUMULATOR_HPP

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace welford {

/**
 * Naive two-pass / sum-of-squares accumulator (used for comparative error benchmarking).
 */
class NaiveAccumulator {
public:
    void update(double x) noexcept;
    [[nodiscard]] double mean() const noexcept;
    [[nodiscard]] double variance() const noexcept;
    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }

private:
    double sum_{0.0};
    double sum_sq_{0.0};
    std::uint64_t count_{0};
};

/**
 * Welford's algorithm for numerically stable online mean and sample variance computation.
 * Avoids catastrophic cancellation inherent in E[X^2] - (E[X])^2.
 */
class WelfordAccumulator {
public:
    void update(double x) noexcept;
    [[nodiscard]] double mean() const noexcept;
    [[nodiscard]] double variance() const noexcept;
    [[nodiscard]] double standard_deviation() const noexcept;
    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }
    void reset() noexcept;

private:
    double mean_{0.0};
    double m2_{0.0};
    std::uint64_t count_{0};
};

/**
 * Online higher-order moment accumulator computing mean, variance, skewness, and excess kurtosis
 * in a single pass up to the 4th central moment.
 */
class MomentsAccumulator {
public:
    void update(double x) noexcept;
    [[nodiscard]] double mean() const noexcept;
    [[nodiscard]] double variance() const noexcept;
    [[nodiscard]] double skewness() const noexcept;
    [[nodiscard]] double excess_kurtosis() const noexcept;
    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }

private:
    double mean_{0.0};
    double m2_{0.0};
    double m3_{0.0};
    double m4_{0.0};
    std::uint64_t count_{0};
};

/**
 * Bivariate online co-variance accumulator for computing running covariance, correlation,
 * and linear regression beta parameter without storing sample history.
 */
class BivariateAccumulator {
public:
    void update(double x, double y) noexcept;
    [[nodiscard]] double covariance() const noexcept;
    [[nodiscard]] double correlation() const noexcept;
    [[nodiscard]] double beta() const noexcept;
    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }

private:
    double mean_x_{0.0};
    double mean_y_{0.0};
    double mxx_{0.0};
    double myy_{0.0};
    double mxy_{0.0};
    std::uint64_t count_{0};
};

/**
 * Exponentially Weighted Moving Average (EWMA) variance estimator.
 */
class EwmaEstimator {
public:
    explicit EwmaEstimator(double decay_lambda) noexcept;
    void update(double x) noexcept;
    [[nodiscard]] double variance() const noexcept;
    [[nodiscard]] double standard_deviation() const noexcept;

private:
    double lambda_{0.94};
    double mean_{0.0};
    double variance_{0.0};
    bool initialized_{false};
};

/**
 * Computes exact reference mean and sample variance using a classic two-pass algorithm.
 */
void compute_reference_stats(const std::vector<double>& data, double& out_mean, double& out_var) noexcept;

/**
 * Generates synthetic benchmark series prone to catastrophic floating-point cancellation.
 */
std::vector<double> generate_cancellation_series(std::size_t count, double base_level);

} // namespace welford

#endif // WELFORD_ACCUMULATOR_HPP
