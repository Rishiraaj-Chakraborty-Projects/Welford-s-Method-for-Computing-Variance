#include "welford/accumulator.hpp"
#include <iostream>
#include <chrono>
#include <random>

int main() {
    constexpr std::size_t NUM_OPERATIONS = 10'000'000;
    std::cout << "Starting Welford Online Accumulator Benchmark (" << NUM_OPERATIONS << " operations)...\n";

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);

    welford::WelfordAccumulator acc;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < NUM_OPERATIONS; ++i) {
        acc.update(dist(rng));
    }
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::nano> duration_ns = end - start;
    const double ns_per_op = duration_ns.count() / static_cast<double>(NUM_OPERATIONS);
    const double ops_per_sec = (static_cast<double>(NUM_OPERATIONS) / duration_ns.count()) * 1e9 / 1e6;

    std::cout << "Results:\n";
    std::cout << "  Total time: " << (duration_ns.count() / 1e6) << " ms\n";
    std::cout << "  Latency:    " << ns_per_op << " ns/op\n";
    std::cout << "  Throughput: " << ops_per_sec << " Million ops/sec\n";
    std::cout << "  Final mean: " << acc.mean() << ", Variance: " << acc.variance() << "\n";

    return 0;
}
