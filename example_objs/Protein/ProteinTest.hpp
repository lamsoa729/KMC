/**********************************************************************
 *                     Unit testing for Protein code                      *
 **********************************************************************/

#include "ProteinBindStatus.hpp"
#include "ProteinConfig.hpp"
#include "ProteinData.hpp"
#include "ProteinType.hpp"

#include "SimToolbox/FDPS/particle_simulator.hpp"
#include "SimToolbox/Sylinder/SylinderNear.hpp"

#include "catch.hpp"

#include <array>
#include <cassert>

TEST_CASE("Test Read ProteinConfig.yaml", "[yamlIO]") {
    ProteinConfig proteinConfig("Protein/ProteinConfig_example.yaml");

    constexpr double small = 1e-6;

    SECTION("Check data read") {
        REQUIRE(proteinConfig.types.size() == 2);
        CHECK(proteinConfig.types[0].fixedEnd0 == true);
        CHECK(abs(proteinConfig.types[0].Ke[0] - 0.246) < small);
        CHECK(abs(proteinConfig.types[1].Ke[0] - 0.246) < small);
        CHECK(abs(proteinConfig.fixedLocations[1][0] - 0.5) < small);
        CHECK(proteinConfig.types[0].LUTablePtr == nullptr);
    }
}
