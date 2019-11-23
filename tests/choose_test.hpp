/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : choose_test
 * @created     : Saturday Nov 23, 2019 11:02:52 MST
 */

#ifndef CHOOSE_TEST_HPP

#define CHOOSE_TEST_HPP
#include "KMC/kmc_choose.hpp"
#include <math.h>

TEST_CASE("Test choosing algorithm based on probability", "[kmc_choose]") {
    // Assemble
    double prob_0, prob_1, no_prob, roll;
    double rolls[3] = {0.0, 0.5, 1.0};
    int activated;

    SECTION("Test no binding probability") {
        prob_0 = prob_1 = 0.;
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == -1);
            CHECK(rolls[i] <= 1.0);
        }
    }
    SECTION("Test 50/50 binding probability with Poisson statistics") {
        prob_0 = prob_1 = .5;
        int heads[3] = {0, 1, -1};
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] <= 1.0);
        }
    }
    SECTION("Test larger than 1 binding probability with Poisson statistics") {
        prob_0 = prob_1 = 1.0;
        int heads[3] = {0, 1, -1};
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] <= 1.0);
        }
    }
    SECTION("Test zeroth head binding probability with Poisson statistics") {
        prob_0 = 1.0;
        prob_1 = 0.0;
        int heads[3] = {0, 0, -1};
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] <= 1.0);
        }
    }
    SECTION("Test first head binding probability with Poisson statistics") {
        prob_0 = 0.0;
        prob_1 = 1.0;
        int heads[3] = {1, 1, -1};
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] <= 1.0);
        }
    }
}

TEST_CASE("Test overflow handling of choose KMC", "[choose_err]") {

    // Assemble
    double prob_0, prob_1, no_prob, roll;
    double rolls[3] = {0.0, 0.5, 1.0};
    double check_rolls[3] = {0.0, 0.5, 1.0};
    int activated;

    SECTION("Overflow values for zeroth head") {
        prob_0 = HUGE_VAL;
        prob_1 = 0.0;
        int heads[3] = {0, 0, 0};
        // Apply
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] == check_rolls[i]); // Rolls should not be changed in
                                               // these cases
        }
    }
    SECTION("Overflow values for first head") {
        prob_0 = 0.0;
        prob_1 = HUGE_VAL;
        int heads[3] = {1, 1, 1};
        // Apply
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] == check_rolls[i]);
        }
    }
    SECTION("Overflow of both heads") {
        prob_0 = HUGE_VAL;
        prob_1 = HUGE_VAL;
        int heads[3] = {0, 1, 1};
        double check_rolls[3] = {0.0, 0.0, 1.0};
        // Apply
        for (int i = 0; i < 3; ++i) {
            // Apply
            activated = choose_kmc_double(prob_0, prob_1, rolls[i]);
            // Assert
            REQUIRE(activated == heads[i]);
            CHECK(rolls[i] == check_rolls[i]);
        }
    }
}
/* TODO: Add in checks for error handling like negative probabilities <23-11-19,
 * ARL> */
// TEST_CASE("Test error handling of choose KMC", "[choose_err]") {

//    // Assemble
//    double prob_0, prob_1, roll;
//    double rolls[3] = {0.0, 0.5, 1.0};
//    double check_rolls[3] = {0.0, 0.5, 1.0};
//    int activated;
//}

#endif /* end of include guard CHOOSE_TEST_HPP */

