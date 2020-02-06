/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : diagnostic_test
 * @created     : Thursday Feb 06, 2020 14:02:37 MST
 */

#ifndef DIAGNOSTIC_TEST_HPP

#define DIAGNOSTIC_TEST_HPP

#include "KMC/two_step_prob.hpp"
#include <cmath>

TEST_CASE("Test two step probability function", "[two_step_prob]") {
    double k1 = 10.;
    double k2 = 20.;
    double t = .1;
    double t1 = .5 * t;

    // Assemble
    SECTION("Symmetric cases") {
        REQUIRE(two_step_prob(k1, k2, t, t1) ==
                Approx(two_step_prob(k2, k1, t, t1)).epsilon(1e-8));
        REQUIRE(two_step_prob(k1, k2, t1) ==
                Approx(two_step_prob(k2, k1, t1)).epsilon(1e-8));
    }
    SECTION("Hard coded cases") {
        REQUIRE(two_step_prob(k1, k1, t) == two_step_prob(k2, k2, t1));
    }
}

TEST_CASE("Test two step probability function", "[two_step_prob]") {
    double k1 = 10.;
    double k2 = 20.;
    double t = .1;
    double t1 = .5 * t;

    // Assemble
    SECTION("Symmetric cases") {
        REQUIRE(two_step_prob(k1, k2, t) ==
                Approx(two_step_prob(k2, k1, t)).epsilon(1e-8));
        REQUIRE(two_step_prob(k1, k2, t1) ==
                Approx(two_step_prob(k2, k1, t1)).epsilon(1e-8));
    }
    SECTION("Hard coded cases") {
        REQUIRE(two_step_prob(k1, k1, t) == two_step_prob(k2, k2, t1));
    }
}

#endif /* end of include guard DIAGNOSTIC_TEST_HPP */

