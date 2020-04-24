/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : diagnostic_test
 * @created     : Thursday Feb 06, 2020 14:02:37 MST
 */

#ifndef DIAGNOSTIC_TEST_HPP

#define DIAGNOSTIC_TEST_HPP

//#include "KMC/two_step_prob.hpp"
#include "KMC/ExampleRod.hpp"
#include "KMC/helpers.hpp"
#include "KMC/integrals.hpp"
#include "KMC/kmc.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/lut_filler_edep.hpp"
#include "KMC/macros.hpp"
#include "catch.hpp"
#include "example_objs/ExampleXlink.hpp"
#include <cmath>

TEST_CASE("Throw error from diagnostic when time step is too large",
          "[diagnostic]") {
    // Assemble
    double KBT = 1.;
    ExampleXlink xlink;
    ExampleRod rod = MockRod<ExampleRod>(0);

    xlink.setMockXlink();

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * rod.radius);
    LookupTable LUT(&lut_filler);

    // int i = 1;
    SECTION("specific case") {
        // Apply
        KMC<ExampleRod> kmc_diag(1., 1., 1., &LUT);
        // Assert
        REQUIRE_THROWS(kmc_diag.Diagnostic(1, 1, 1, 1));
    }
}

#endif /* end of include guard DIAGNOSTIC_TEST_HPP */

