# Welford online moments

Single-pass, numerically stable mean, variance, skewness, kurtosis, covariance and exponentially weighted variance. Written in C++17, no dependencies.

The point of the repo is one number.

```
100000 observations of size 1e+08 with a spread of 1e-3 around it

  method                               variance      rel error
  two pass (reference)         4.0000062099e-06              -
  naive sum of squares         3.8320000000e+03      9.580e+08
  Welford                      4.0000063225e-06      2.816e-08
```

The naive estimator is wrong by a factor of a billion. Not noisy. Wrong.

## Why it happens

Everyone learns `Var(X) = E[X^2] - E[X]^2` and it is correct in exact arithmetic. In floating point it computes the difference of two large, nearly equal quantities.

Here the data sit at `1e8` with a spread of `1e-3`. So `sumsq/n` is about `1e16` and `(sum/n)^2` is about `1e16`, and they agree to roughly sixteen significant figures. A double carries about sixteen. The answer lives entirely in the digits that were rounded away, so the subtraction returns rounding noise.

Welford never forms either large quantity. It carries the running mean and the running sum of squared deviations, updating both with

```
n    += 1
delta = x - mean
mean += delta / n
m2   += delta * (x - mean)
```

Note that the second `mean` in the last line is the updated one. That asymmetry is not a typo and swapping it breaks the algorithm.

## Why this matters outside a textbook

Financial data has exactly this shape. An index level of 5,000 with daily moves of 0.5 points. A bond priced at 99.87 moving in basis points. A futures contract at 72,000 with a tick of 5. Any time the signal is small relative to the level, the naive formula degrades, and it degrades silently. It does not throw, it does not warn, it returns a plausible-looking number.

The failure is worst in exactly the situation you care about: low volatility. High volatility data has a large spread relative to the level and the naive formula survives it. Calm markets are where it breaks.

## What else is in here

**Fourth-moment extension.** Skewness and kurtosis from the same single pass, using the standard Pébay update. Excess kurtosis matters for options work because it is the number that says how far the return distribution sits from lognormal.

**Online covariance.** Two-variable Welford, which gives you covariance, correlation and the regression slope in one pass. That slope is the minimum-variance hedge ratio, so this is the streaming version of a beta calculation.

**Exponentially weighted variance.** What risk systems actually run, because the equally weighted estimator treats a return from two years ago the same as yesterday's. The demo puts a regime break in the middle of the sample and shows the tradeoff: faster decay tracks the break sooner and is noisier afterward.

## Build

PowerShell:

```powershell
g++ -std=c++17 -O2 -o welford.exe src\main.cpp
.\welford.exe
```

CMake, any platform:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Make:

```
make && make run
```

## Run

```
welford                      full demonstration
welford --test               self checks, exit code 0 or 1
welford --n 1000000          more observations
welford --level 1e12         make the cancellation worse
welford --help
```

Push `--level` up and watch the naive error grow. At `1e12` there is nothing left at all.

## Self test

Ten checks, wired into CTest so CI runs them.

```
  mean of the textbook sample                    ok
  sample variance of the same                    ok
  Welford within 1e-6 relative on hard data      ok
  naive blows up by more than 1e6 relative       ok
  uniform skewness near zero                     ok
  uniform excess kurtosis near -1.2              ok
  perfect linear pair gives beta 3               ok
  perfect linear pair gives corr 1               ok
  variance is shift invariant                    ok
  single observation reports zero variance       ok
```

The shift-invariance check is the one that matters. Variance does not change when you add a constant to every observation, so an estimator whose answer moves when you shift the data is broken. The naive one fails that test badly and Welford passes it to nine digits.

The uniform kurtosis check is a real test rather than a tautology: a uniform distribution has excess kurtosis of exactly `-1.2`, so recovering it confirms the fourth-moment update is right and not merely self-consistent.

## Notes

`Naive` is kept in the source deliberately. It is not dead code, it is the control.

The variance returned uses the `n-1` denominator. It is unbiased for the variance but the square root of it is still biased low for the standard deviation, by roughly `sigma/(4n)`, since the square root is concave. Small, but real if you are estimating volatility from twenty observations.

## Licence

MIT.
