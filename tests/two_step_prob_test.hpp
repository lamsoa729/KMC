/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : diagnostic_test
 * @created     : Thursday Feb 06, 2020 14:02:37 MST
 */

#ifndef TWO_STEP_PROB_TEST_HPP

#define TWO_STEP_PROB_TEST_HPP

#include "KMC/two_step_prob.hpp"
#include <cmath>

TEST_CASE("Test time at which max probability is reached between two steps",
          "[time_max]") {
    // Assemble
    double k1 = 10.;
    double k2 = 20.;
    double t = .1;
    double t1 = .5 * t;

    SECTION("Hard coded cases") {
        REQUIRE(max_time_of_two_step_prob(k1, k2, t) ==
                Approx(0.05491770806462643).epsilon(1e-8));
        REQUIRE(max_time_of_two_step_prob(k2, k1, t) ==
                Approx(t - 0.05491770806462643).epsilon(1e-8));
        REQUIRE(max_time_of_two_step_prob(10 * k1, k2, t) ==
                Approx(0.02831502220553105).epsilon(1e-8));
        REQUIRE(max_time_of_two_step_prob(k2, 10 * k1, t) ==
                Approx(t - 0.02831502220553105).epsilon(1e-8));
    }
    SECTION("Always true results of max time") {
        REQUIRE(max_time_of_two_step_prob(k1, k1, t) ==
                Approx(t1).epsilon(1e-8));
        REQUIRE(max_time_of_two_step_prob(k2, k2, t) ==
                Approx(t1).epsilon(1e-8));
    }
}

TEST_CASE("Test two step probability function", "[two_step_prob]") {
    double k1 = 10.;
    double k2 = 20.;
    double t = .1;
    double t1 = .75 * t;

    // Assemble
    SECTION("Symmetry cases") {
        REQUIRE(two_step_prob(k1, k2, t, t1) ==
                Approx(two_step_prob(k2, k1, t, t - t1)).epsilon(1e-8));
    }

    SECTION("Zero cases") {
        REQUIRE(two_step_prob(0., k1, t, t1) == 0.);
        REQUIRE(two_step_prob(k2, 0., t, t1) == 0.);
        REQUIRE(two_step_prob(k2, k1, 0., 0.) == 0.);
    }

    SECTION("Scaling cases") {
        REQUIRE(two_step_prob(k1, k1, t, t1) ==
                Approx(two_step_prob(k2, k2, .5 * t, .5 * t1)).epsilon(1e-8));
    }
}

TEST_CASE("Test two step max probability function", "[two_step_max_prob]") {
    double k1 = 10.;
    double k2 = 20.;
    double t = .1;

    // Assemble
    SECTION("Symmetry cases") {
        REQUIRE(two_step_max_prob(k1, k2, t) ==
                Approx(two_step_max_prob(k2, k1, t)).epsilon(1e-8));
        REQUIRE(two_step_max_prob(10 * k1, k2, t) ==
                Approx(two_step_max_prob(k2, 10 * k1, t)).epsilon(1e-8));
        REQUIRE(two_step_max_prob(k1, k2, .5 * t) ==
                Approx(two_step_max_prob(k2, k1, .5 * t)).epsilon(1e-8));
        REQUIRE(two_step_max_prob(k1, k2, 10 * t) ==
                Approx(two_step_max_prob(k2, k1, 10 * t)).epsilon(1e-8));
    }
}

#endif /* end of include guard TWO_STEP_PROB_TEST_HPP */

