# Online Variance & High-Order Moment Estimation (Welford's Algorithm)

A C++17 library implementing numerically stable online algorithms for streaming variance, higher-order central moments (skewness and excess kurtosis), bivariate covariance/correlation, and exponentially weighted moving average (EWMA) estimations in a single pass.

## Overview

In quantitative trading applications and high-frequency streaming telemetry, computing sample variance using the standard textbook formula ($E[X^2] - (E[X])^2$) fails when processing large inputs or streams with small variance around a high mean. Floating-point cancellation destroys precision in the subtraction of two large floating-point numbers.

Welford's algorithm computes running mean and sample variance incrementally, avoiding catastrophic cancellation and requiring $O(1)$ memory space without retaining sample history.

## Mathematics & Algorithm Specification

### Single Pass Mean & Variance
Given sample stream $x_1, x_2, \dots, x_n$:

$$\bar{x}_n = \bar{x}_{n-1} + \frac{x_n - \bar{x}_{n-1}}{n}$$

$$M_{2,n} = M_{2,n-1} + (x_n - \bar{x}_{n-1})(x_n - \bar{x}_n)$$

$$s^2 = \frac{M_{2,n}}{n - 1} \quad (n > 1)$$

### Higher-Order Moments (Skewness & Kurtosis)
Extending the recurrence up to the 4th central moment allows tracking skewness $g_1$ and excess kurtosis $g_2$ online in single pass without store-and-revisit.

### Bivariate Covariance & Regression Beta
Tracks running covariance between two streaming signals $(x_i, y_i)$:

$$C_n = C_{n-1} + (x_n - \bar{x}_n)(y_n - \bar{y}_{n-1})$$

$$\text{Cov}(X, Y) = \frac{C_n}{n - 1}, \quad \beta = \frac{C_n}{M_{2,X}}$$

## Repository Structure

```text
01-welford-online-variance/
├── CMakeLists.txt
├── README.md
├── include/
│   └── welford/
│       └── accumulator.hpp    # Class declarations & API
├── src/
│   └── accumulator.cpp        # Welford accumulation implementation
├── tests/
│   └── test_welford.cpp       # Numerical stability unit test suite
├── benchmarks/
│   └── bench_welford.cpp      # Microbenchmark runner
└── apps/
    └── main.cpp               # Interactive CLI demo application
```

## Building & Testing

### Prerequisites
* C++17 compatible compiler (`g++`, `clang++`, or MSVC)
* CMake 3.15+

### Build Instructions

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Running Executables

```bash
# Run unit test suite
./welford_tests

# Run benchmark harness
./welford_benchmark

# Run CLI demonstration tool
./welford_demo --n 100000 --level 1e8
```

## Numerical Precision Comparison

| Method | Sample Variance ($s^2$) | Relative Error |
| :--- | :--- | :--- |
| Two-Pass Reference | `2.9166666667e-07` | Baseline (`0.0`) |
| Naive Sum of Squares | Imprecise / Loss of precision | $> 10^6$ (Loss of Significance) |
| **Welford (One-Pass)** | `2.9166666667e-07` | $< 10^{-14}$ |

## License

MIT License.
