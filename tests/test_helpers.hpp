/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : test_helpers
 * @created     : Thursday Jan 30, 2020 16:44:22 MST
 */

#ifndef TEST_HELPERS_HPP

#define TEST_HELPERS_HPP

#include <cmath>
constexpr double ABSTOL = 1e-3;
constexpr double RELTOL = 1e-2;

constexpr double REVERSEFAC = 10;
/**
 * EXPLAIN:
 * Reverse lookup on regions where the CDF is flat gives large error
 * But this is a rare event in simulations because this requires
 * the U01 rng to be very close to either 0 or 1
 * WARNING:
 * When spring is stiff and distPerp is large,
 * the reverse lookup error is very large
 * because function values are tiny and very flat.
 */

inline double absError(double a, double b) { return fabs(a - b); }

inline double relError(double a, double b) {
    return b < 1e-8 ? absError(a, b) : fabs((a - b) / b);
}

inline bool errorPass(double a, double b, double fac = 1) {
    return absError(a, b) < fac * ABSTOL || relError(a, b) < fac * RELTOL;
}

#endif /* end of include guard TEST_HELPERS_HPP */

