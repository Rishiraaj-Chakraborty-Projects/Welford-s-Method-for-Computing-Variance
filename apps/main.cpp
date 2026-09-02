#include "welford/accumulator.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

void print_usage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [--n N] [--level L]\n";
    std::cout << "  --n N      Number of observation samples (default: 100000)\n";
    std::cout << "  --level L  Base offset level for cancellation testing (default: 1e8)\n";
}

int main(int argc, char** argv) {
    std::size_t num_samples = 100000;
    double level = 1e8;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--n" && i + 1 < argc) {
            num_samples = static_cast<std::size_t>(std::stoul(argv[++i]));
        } else if (arg == "--level" && i + 1 < argc) {
            level = std::stod(argv[++i]);
        }
    }

    std::cout << "========================================================\n";
    std::cout << "   WELFORD ONLINE ACCUMULATOR DEMONSTRATION SUITE       \n";
    std::cout << "========================================================\n\n";

    const auto series = welford::generate_cancellation_series(num_samples, level);

    welford::NaiveAccumulator naive_acc;
    welford::WelfordAccumulator welford_acc;
    welford::MomentsAccumulator moments_acc;

    for (double val : series) {
        naive_acc.update(val);
        welford_acc.update(val);
        moments_acc.update(val);
    }

    double ref_mean = 0.0, ref_var = 0.0;
    welford::compute_reference_stats(series, ref_mean, ref_var);

    std::cout << "Series Specification: " << num_samples << " observations around mean level " << level << "\n\n";

    std::cout << std::left << std::setw(28) << "Method"
              << std::right << std::setw(22) << "Variance"
              << std::setw(18) << "Relative Error" << "\n";
    std::cout << std::string(68, '-') << "\n";

    std::cout << std::left << std::setw(28) << "Reference (Two-Pass)"
              << std::right << std::setw(22) << std::scientific << std::setprecision(10) << ref_var
              << std::setw(18) << "-" << "\n";

    const double naive_rel_err = std::abs(naive_acc.variance() - ref_var) / ref_var;
    std::cout << std::left << std::setw(28) << "Naive (Sum of Squares)"
              << std::right << std::setw(22) << std::scientific << std::setprecision(10) << naive_acc.variance()
              << std::setw(18) << std::scientific << std::setprecision(3) << naive_rel_err << "\n";

    const double welford_rel_err = std::abs(welford_acc.variance() - ref_var) / ref_var;
    std::cout << std::left << std::setw(28) << "Welford (One-Pass)"
              << std::right << std::setw(22) << std::scientific << std::setprecision(10) << welford_acc.variance()
              << std::setw(18) << std::scientific << std::setprecision(3) << welford_rel_err << "\n";

    const double moments_rel_err = std::abs(moments_acc.variance() - ref_var) / ref_var;
    std::cout << std::left << std::setw(28) << "Moments (4th Order)"
              << std::right << std::setw(22) << std::scientific << std::setprecision(10) << moments_acc.variance()
              << std::setw(18) << std::scientific << std::setprecision(3) << moments_rel_err << "\n";

    std::cout << "\nHigher Central Moments (Single Pass):\n";
    std::cout << "  Mean:            " << std::fixed << std::setprecision(6) << moments_acc.mean() << "\n";
    std::cout << "  Sample Variance: " << std::scientific << std::setprecision(6) << moments_acc.variance() << "\n";
    std::cout << "  Skewness:        " << std::fixed << std::setprecision(6) << moments_acc.skewness() << "\n";
    std::cout << "  Excess Kurtosis: " << std::fixed << std::setprecision(6) << moments_acc.excess_kurtosis() << "\n";

    return 0;
}
