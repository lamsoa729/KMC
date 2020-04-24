/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : diagnostic_test
 * @created     : Thursday Feb 06, 2020 14:02:37 MST
 */

#ifndef DIAGNOSTIC_TEST_HPP

#define DIAGNOSTIC_TEST_HPP

#include "KMC/ExampleRod.hpp"
#include "KMC/helpers.hpp"
#include "KMC/kmc.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/lut_filler_edep.hpp"
#include "KMC/macros.hpp"
#include "catch.hpp"
#include <cmath>

TEST_CASE("Throw error from non-energy(force) dependent unbinding diagnostic",
          "[ne-diagnostic]") {
    // Assemble
    double dt = 1.;
    double diffConst = 0; // Small because we want to use r_cutoff = .5
    double r_cutoff = .5;
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(M_PI, 0, 1); // Gaussian distribution normalized to 1
    LookupTable LUT(&lut_filler);

    // Apply
    KMC<ExampleRod> kmc_diag(r_cutoff, diffConst, dt, &LUT);

    // Assert
    REQUIRE_NOTHROW(kmc_diag.Diagnostic(.0001, .0001, .0001, .0001));
    SECTION("U->S->U testing") {
        // Assert
        REQUIRE_THROWS(kmc_diag.Diagnostic(1, 1, .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1000., 1., .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1., 1000., .0001, .0001));
    }
    SECTION("U->S->D testing") {
        // Assert
        REQUIRE_THROWS(kmc_diag.Diagnostic(1, 1, .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1000., 1., .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1., 1000., .0001, .0001));
    }
    SECTION("D->S->D testing") {
        // Assert
        REQUIRE_THROWS(kmc_diag.Diagnostic(1, 1, .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1000., 1., .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1., 1000., .0001, .0001));
    }
    SECTION("D->S->U testing") {
        // Assert
        REQUIRE_THROWS(kmc_diag.Diagnostic(1, 1, .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1000., 1., .0001, .0001));
        REQUIRE_THROWS(kmc_diag.Diagnostic(1., 1000., .0001, .0001));
    }
}

TEST_CASE("Throw error from energy(force) dependent unbinding diagnostic",
          "[unbind-dep-diagnostic]") {
    // Assemble
    double dt = 1.;
    double diffConst = 0;
    double r_cutoff = .5;
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(M_PI, 0, 1); // Gaussian distribution normalized to 1
    LookupTable LUT(&lut_filler);

    SECTION("U->S->U testing") {
        // Apply
        KMC<ExampleRod> kmc_diag(r_cutoff, diffConst, dt, &LUT);
        // Assert
        REQUIRE_NOTHROW(kmc_diag.DiagnosticUnBindDep(.0001, .0001, .0001));
        REQUIRE_THROWS(kmc_diag.DiagnosticUnBindDep(1, 1, .0001));
        REQUIRE_THROWS(kmc_diag.DiagnosticUnBindDep(1000., 1., .0001));
        REQUIRE_THROWS(kmc_diag.DiagnosticUnBindDep(1., 1000., .0001));
    }
}

#endif /* end of include guard DIAGNOSTIC_TEST_HPP */

