#include "welford/accumulator.hpp"
#include <cmath>
#include <algorithm>

namespace welford {

// --- NaiveAccumulator ---
void NaiveAccumulator::update(double x) noexcept {
    sum_ += x;
    sum_sq_ += x * x;
    ++count_;
}

double NaiveAccumulator::mean() const noexcept {
    return count_ > 0 ? sum_ / static_cast<double>(count_) : 0.0;
}

double NaiveAccumulator::variance() const noexcept {
    if (count_ == 0) return 0.0;
    const double m = mean();
    return (sum_sq_ / static_cast<double>(count_)) - (m * m);
}

// --- WelfordAccumulator ---
void WelfordAccumulator::update(double x) noexcept {
    ++count_;
    const double delta = x - mean_;
    mean_ += delta / static_cast<double>(count_);
    m2_ += delta * (x - mean_);
}

double WelfordAccumulator::mean() const noexcept {
    return mean_;
}

double WelfordAccumulator::variance() const noexcept {
    return count_ > 1 ? m2_ / static_cast<double>(count_ - 1) : 0.0;
}

double WelfordAccumulator::standard_deviation() const noexcept {
    return std::sqrt(variance());
}

void WelfordAccumulator::reset() noexcept {
    mean_ = 0.0;
    m2_ = 0.0;
    count_ = 0;
}

// --- MomentsAccumulator ---
void MomentsAccumulator::update(double x) noexcept {
    const double n1 = static_cast<double>(count_);
    ++count_;
    const double n = static_cast<double>(count_);
    
    const double delta = x - mean_;
    const double delta_n = delta / n;
    const double delta_n2 = delta_n * delta_n;
    const double term1 = delta * delta_n * n1;

    mean_ += delta_n;
    m4_ += term1 * delta_n2 * (n * n - 3.0 * n + 3.0) + 6.0 * delta_n2 * m2_ - 4.0 * delta_n * m3_;
    m3_ += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * m2_;
    m2_ += term1;
}

double MomentsAccumulator::mean() const noexcept {
    return mean_;
}

double MomentsAccumulator::variance() const noexcept {
    return count_ > 1 ? m2_ / static_cast<double>(count_ - 1) : 0.0;
}

double MomentsAccumulator::skewness() const noexcept {
    if (count_ <= 2 || m2_ <= 0.0) return 0.0;
    return std::sqrt(static_cast<double>(count_)) * m3_ / std::pow(m2_, 1.5);
}

double MomentsAccumulator::excess_kurtosis() const noexcept {
    if (count_ <= 3 || m2_ <= 0.0) return 0.0;
    return (static_cast<double>(count_) * m4_) / (m2_ * m2_) - 3.0;
}

// --- BivariateAccumulator ---
void BivariateAccumulator::update(double x, double y) noexcept {
    ++count_;
    const double n = static_cast<double>(count_);
    const double delta_x = x - mean_x_;
    const double delta_y = y - mean_y_;

    mean_x_ += delta_x / n;
    mean_y_ += delta_y / n;
    mxx_ += delta_x * (x - mean_x_);
    myy_ += delta_y * (y - mean_y_);
    mxy_ += delta_x * (y - mean_y_);
}

double BivariateAccumulator::covariance() const noexcept {
    return count_ > 1 ? mxy_ / static_cast<double>(count_ - 1) : 0.0;
}

double BivariateAccumulator::correlation() const noexcept {
    const double denom = std::sqrt(mxx_ * myy_);
    return denom > 0.0 ? mxy_ / denom : 0.0;
}

double BivariateAccumulator::beta() const noexcept {
    return mxx_ > 0.0 ? mxy_ / mxx_ : 0.0;
}

// --- EwmaEstimator ---
EwmaEstimator::EwmaEstimator(double decay_lambda) noexcept
    : lambda_(decay_lambda) {}

void EwmaEstimator::update(double x) noexcept {
    if (!initialized_) {
        mean_ = x;
        variance_ = 0.0;
        initialized_ = true;
        return;
    }
    const double delta = x - mean_;
    mean_ += (1.0 - lambda_) * delta;
    variance_ = lambda_ * variance_ + (1.0 - lambda_) * delta * delta;
}

double EwmaEstimator::variance() const noexcept {
    return variance_;
}

double EwmaEstimator::standard_deviation() const noexcept {
    return std::sqrt(variance_);
}

// --- Helper Functions ---
void compute_reference_stats(const std::vector<double>& data, double& out_mean, double& out_var) noexcept {
    if (data.empty()) {
        out_mean = 0.0;
        out_var = 0.0;
        return;
    }
    double sum = 0.0;
    for (double val : data) sum += val;
    out_mean = sum / static_cast<double>(data.size());

    double sum_sq_diff = 0.0;
    for (double val : data) {
        const double diff = val - out_mean;
        sum_sq_diff += diff * diff;
    }
    out_var = data.size() > 1 ? sum_sq_diff / static_cast<double>(data.size() - 1) : 0.0;
}

std::vector<double> generate_cancellation_series(std::size_t count, double base_level) {
    std::vector<double> series;
    series.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        series.push_back(base_level + (static_cast<double>(i % 7) - 3.0) * 1e-3);
    }
    return series;
}

} // namespace welford
