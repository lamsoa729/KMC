/**********************************************************************
 *                     Unit testing for KMC<ExampleRod> code *
 **********************************************************************/

#include "KMC/ExampleRod.hpp"
#include "KMC/helpers.hpp"
#include "KMC/kmc.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/macros.hpp"
#include "example_objs/ExampleXlink.hpp"

#include "catch.hpp"

#include <array>
#include <cassert>

template <class TRod>
TRod MockRod(int id) {
    TRod rod;
    rod.gid = id;
    rod.length = 400;
    rod.rank = 0;
    rod.radius = .5;
    for (int i = 0; i < 3; ++i) {
        rod.pos[i] = 0;
        rod.direction[i] = 0;
    }
    rod.direction[0] = 1;
    return rod;
}

TEST_CASE("Test CalcProbUS for KMC<ExampleRod> class", "[calc01prob]") {
    // Time step
    double dt = .0001;
    double KBT = 1.;
    // Protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // Sylinder data
    ExampleRod rod = MockRod<ExampleRod>(0);

    SECTION("Test binding probability when crosslinker is in center of rod") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0);
        REQUIRE(ABS(prob - 0.0001909859317102744) <= 1e-8);
    }
    SECTION(
        "Test binding probability when crosslinker is above and outside rod ") {
        xlink.pos[2] = 2.0;
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 2.0);
        CHECK(prob == 0);
    }
    SECTION(
        "Test binding probability when crosslinker is at plus end of rod.") {
        double pos[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
            pos[i] = rod.pos[i] + (rod.direction[i] * rod.length * .5);
        }
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == rod.length * .5);
        CHECK(kmc.getDistMin(0) == 0);
        CHECK(ABS(prob - 0.0000954929658551372) <= 1e-8);
    }
    SECTION(
        "Test binding probability when crosslinker is at minus end of rod.") {
        double pos[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
            pos[i] = rod.pos[i] - (rod.direction[i] * rod.length * .5);
        }
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == rod.length * -.5);
        CHECK(kmc.getDistMin(0) == 0);
        REQUIRE(ABS(prob - 0.0000954929658551372) <= 1e-8);
    }
    SECTION("Test binding probability when crosslinker is on surface of rod "
            "above its center.") {
        double pos[3] = {0, 0, 0.5};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.5);
        REQUIRE(prob == 0);
    }
    SECTION("Test binding probability when crosslinker is on surface of rod's "
            "tip.") {
        double pos[3] = {rod.length * 0.5, 0, 0.5};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == rod.length * .5);
        CHECK(kmc.getDistMin(0) == 0.5);
        REQUIRE(prob == 0);
    }
    SECTION("Test binding probability when crosslinker is between surface of "
            "rod  and rod axis above its center.") {
        double pos[3] = {0, 0, 0.25};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.25);
        REQUIRE(ABS(prob - 0.0001653986686265376) <= 1e-8);
    }
    SECTION("Test binding probability when crosslinker binds to a rod smaller "
            "than cutoff length") {
        double pos[3] = {0, 0, 0};
        rod.length = .5;
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
        double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.0);
        REQUIRE(ABS(prob - 0.0000954929658551372) < SMALL);
    }
}

TEST_CASE("Test CalcProbSU for KMC class", "[calc_prob_su]") {
    // Time step
    double dt = .0001;
    double KBT = 1.;
    // Protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // Sylinder data
    ExampleRod rod = MockRod<ExampleRod>(1);

    SECTION("Test binding is working properly") {
        int end_bound = 0;
        double end_loc = 0;
        xlink.setBind(end_bound, rod.gid, rod.direction, rod.pos, end_loc,
                      rod.length, rod.rank);
        xlink.updateGeometryWithBind();

        CHECK(xlink.idBind[0] == rod.gid);
        CHECK(xlink.rankBind[0] == rod.rank);
        CHECK(xlink.lenBind[0] == rod.length);
        CHECK(xlink.directionBind[0][0] == rod.direction[0]);
        CHECK(xlink.directionBind[0][1] == rod.direction[1]);
        CHECK(xlink.directionBind[0][2] == rod.direction[2]);
        CHECK(xlink.posEndBind[0][0] == rod.pos[0]);
        CHECK(xlink.posEndBind[0][1] == rod.pos[1]);
        CHECK(xlink.posEndBind[0][2] == rod.pos[2]);
        CHECK(xlink.centerBind[0][0] == rod.pos[0]);
        CHECK(xlink.centerBind[0][1] == rod.pos[1]);
        CHECK(xlink.centerBind[0][2] == rod.pos[2]);
        CHECK(xlink.pos[0] == rod.pos[0]);
        CHECK(xlink.pos[1] == rod.pos[1]);
        CHECK(xlink.pos[2] == rod.pos[2]);
        CHECK(xlink.idBind[1] == IDUB);
    }

    SECTION("Test unbinding probability") {
        xlink.setBind(0, rod.gid, rod.direction, rod.pos, 0, rod.length,
                      rod.rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 1, xlink.getRcutUS());
        kmc.CalcProbSU(xlink.getUnbindingFactorDS(0, dt, KBT));

        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
    }
}

TEST_CASE("Test CalcProbSD for KMC<ExampleRod> class", "[calc_prob_sd]") {
    double dt = .0001;
    double KBT = 1.;
    // Protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // //ProteinType *pType = &xlink.property;
    // Sylinder data
    ExampleRod rod0, rod1;
    rod0 = MockRod<ExampleRod>(0);
    rod0.length = 20;
    rod1 = MockRod<ExampleRod>(1);
    rod1.length = 20;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * rod0.radius);
    xlink.LUTablePtr = &LUT;

    SECTION("Test binding to 2 parallel overlapping rods") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0004899204301239162) < SMALL);
    }
    SECTION("Test binding to 2 parallel vertically separated rods.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Change the position of rod1 vertically
        rod1.pos[0] = 0;
        rod1.pos[1] = 1.0;
        rod1.pos[2] = 0;
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0005481166497130916) < SMALL);
    }
    SECTION("Test binding to 2 parallel rods with crosslinker bound at rod "
            "plus end.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Change the position of rod1 vertically
        rod1.pos[0] = 0;
        rod1.pos[1] = 1.0;
        rod1.pos[2] = 0;
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 1, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0002740583248565458) < SMALL);
    }
    SECTION("Test binding to 2 parallel rods with crosslinker bound at rod "
            "minus end.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = -rod0.length * .5;
        // Change the position of rod1 vertically
        rod1.pos[0] = 0;
        rod1.pos[1] = 1.0;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        std::cout << xlink.getRcutSD() << std::endl;
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0002740583248565458) < SMALL);
    }
    SECTION("Test binding to 2 parallel inline rods with crosslinker bound at "
            "rod plus end and other rod shifted by rod length + D.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        // Change the position of rod1 horizontally
        rod1.pos[0] = 21.0;
        rod1.pos[1] = 0;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0002108938529207641) < SMALL);
    }
    SECTION("Test binding to 2 parallel with crosslinker bound at "
            "rod plus end and other rod shifted horizontially by rod length + "
            "D and vertically by D") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        // Change the position of rod1 horizontally
        rod1.pos[0] = 21.0;
        rod1.pos[1] = 1.0;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.000204684512149819) < SMALL);
    }
    SECTION("Test binding to 2 perpendicular rods crosslinked in the center "
            "and one rod shifted over by D.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        // Change the position and direction of rod0
        rod0.pos[0] = 0;
        rod0.pos[1] = 1.0;
        rod0.pos[2] = 0;
        rod0.direction[0] = 0;
        rod0.direction[1] = 1.0;
        rod0.direction[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                           xlink.freeLength, xlink.getBindingFactorSD(1, dt));
        // Check to make sure that probability matches up
        CHECK(ABS(prob - 0.0005481166497130916) < SMALL);
    }
}

TEST_CASE("Test CalcProbDS for KMC<ExampleRod> class", "[calc21prob]") {
    double dt = .0001;
    double KBT = 1.;
    // Protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    // Sylinder data
    ExampleRod rod0, rod1;
    rod0 = MockRod<ExampleRod>(0);
    rod0.length = 20;
    rod1 = MockRod<ExampleRod>(1);
    rod1.length = 20;
    rod1.pos[0] = 0;
    rod1.pos[1] = 1.0;
    rod1.pos[2] = 0;

    // Set up lookup table
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * rod0.radius);
    xlink.LUTablePtr = &LUT;

    SECTION("Test unbinding of end 0 from rod0.") {
        double end_loc = 0;
        xlink.setBind(0, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        xlink.setBind(1, rod1.gid, rod1.direction, rod1.pos, end_loc,
                      rod0.length, rod0.rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(0, dt, KBT));
        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
    }
    SECTION("Test unbinding of end 1 from rod1.") {
        double end_loc = 0;
        xlink.setBind(0, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        xlink.setBind(1, rod1.gid, rod1.direction, rod1.pos, end_loc,
                      rod0.length, rod0.rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 2, xlink.getRcutSD());
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, dt, KBT));
        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
    }
}

TEST_CASE("Test binding probabilities of one head when near multiple rods.",
          "[calc01tot_probs]") {
    double dt = .0001;
    double KBT = 1.;
    ExampleRod *ep_j[4];
    for (int i = 0; i < 4; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
    }
    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * ep_j[0]->radius);
    xlink.LUTablePtr = &LUT;
    std::vector<int> Uniquefilter{1, 1, 1, 1};
    SECTION("Two semi-overlapping parallel rods with protein in between.") {
        ep_j[0]->pos[0] = 0.0;
        ep_j[0]->pos[1] = .5;
        ep_j[0]->pos[2] = 0.0;
        double pos[3] = {0, 0.25, 0};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 2, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        REQUIRE(ABS(kmc.getTotProb() - 0.0003307973372530752) < SMALL);
    }
    SECTION("Test to make sure filter is working.") {
        // Set filters for rods 2 and 3 to 0
        Uniquefilter[2] = Uniquefilter[3] = 0;
        ep_j[0]->pos[0] = 0.0;
        ep_j[0]->pos[1] = 0.5;
        ep_j[0]->pos[2] = 0.0;
        double pos[3] = {0, 0.25, 0};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 4, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        REQUIRE(ABS(kmc.getTotProb() - 0.0003307973372530752) < SMALL);
    }
    SECTION("4 parallel rods with protein in the center.") {
        ep_j[1]->pos[0] = 0.0;
        ep_j[1]->pos[1] = .5;
        ep_j[1]->pos[2] = 0.0;

        ep_j[2]->pos[0] = 0.0;
        ep_j[2]->pos[1] = 0.25;
        ep_j[2]->pos[2] = 0.25;

        ep_j[3]->pos[0] = 0.0;
        ep_j[3]->pos[1] = 0.25;
        ep_j[3]->pos[2] = -0.25;

        double pos[3] = {0, 0.25, 0};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 4, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        CHECK(ABS(kmc.getTotProb() - (0.0006615946745061504)) < SMALL);
    }
    // Clean up testing
    for (int i = 0; i < 4; ++i) {
        delete ep_j[i];
    }
}

TEST_CASE("Test binding probabilities of one head when near multiple rods and "
          "other head is bound.",
          "[calc12tot_probs]") {
    double dt = .0001;
    double KBT = 1.;
    ExampleRod *ep_j[4];
    for (int i = 0; i < 4; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
    }
    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * ep_j[0]->radius);
    xlink.LUTablePtr = &LUT;

    std::vector<int> Uniquefilter{1, 1, 1, 1};

    SECTION("No rods nearby") {
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 1, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb()) < SMALL);
    }
    SECTION("One rod nearby and above bound crosslinker") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.0005481166497130916) < SMALL);
    }
    SECTION("Two rods nearby, above and below bound crosslinker") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.001096233299426183) < SMALL);
    }
    SECTION("Two rods nearby, above and perpendicular but same distance from "
            "bound crosslinker.") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 1.0;
        ep_j[2]->pos[1] = 1.0;
        ep_j[2]->pos[2] = 0;
        ep_j[2]->direction[0] = 0;
        ep_j[2]->direction[1] = 1;
        ep_j[2]->direction[2] = 0;

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.001096233299426183) < SMALL);
    }
    SECTION("Two rods nearby, above and forward touching cutoff.") {
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = xlink.posEndBind[0][1] + xlink.getRcutSD();
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = ep_j[0]->length + xlink.getRcutSD();
        ep_j[2]->pos[1] = 0;
        ep_j[2]->pos[2] = 0;
        // ep_j[2]->setpos(pos);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb()) <
              dt * .001); // 3 decimal accuracy more than dt
    }
    SECTION("Two rods nearby, above and below outside cutoff") {
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);

        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = xlink.posEndBind[0][1] + xlink.getRcutSD() + 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = xlink.posEndBind[0][1] - xlink.getRcutSD() - 1.0;
        ep_j[2]->pos[2] = 0;

        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));

        CHECK(ABS(kmc.getTotProb()) < SMALL);
    }
    SECTION("Three rods nearby, above and below bound crosslinker. Two "
            "overlapping.") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;

        ep_j[3]->pos[0] = 0;
        ep_j[3]->pos[1] = -1.0;
        ep_j[3]->pos[2] = 0;
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.001644349949139275) < SMALL);
    }
    SECTION("Three rods nearby, above and below bound crosslinker. Two "
            "overlapping but one filtered.") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;

        ep_j[3]->pos[0] = 0;
        ep_j[3]->pos[1] = -1.0;
        ep_j[3]->pos[2] = 0;
        Uniquefilter[3] = 0;
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.001096233299426183) < SMALL);
    }
    SECTION("Three rods nearby above, below, and on, all inside rcut.") {
        ep_j[1]->length = ep_j[2]->length = ep_j[3]->length = 2.0;

        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.0003456279751132963) < SMALL);
    }
    SECTION("Three rods partially inside rcut with different orientations.") {
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        ep_j[1]->pos[0] = 0.5 * ep_j[1]->length;
        ep_j[1]->pos[1] = 0.0;
        ep_j[1]->pos[2] = 0.0;

        ep_j[2]->pos[0] = -0.5 * ep_j[2]->length;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0.0;
        // pos[0] = 0;
        ep_j[3]->direction[0] = 0.0;
        ep_j[3]->direction[1] = 1.0;
        ep_j[3]->direction[2] = 0.0;
        // pos[1] = 0.5 * ep_j[3]->length + 1.0;
        // ep_j[3]->setpos(pos);
        ep_j[3]->pos[0] = 0;
        ep_j[3]->pos[1] = 0.5 * ep_j[3]->length + 1.0;
        ep_j[3]->pos[2] = 0.0;
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD());

        kmc.CalcTotProbsSD(ep_j, Uniquefilter, ep_j[0]->gid, xlink.lambda,
                           xlink.kappa, 1. / KBT, xlink.freeLength,
                           xlink.getBindingFactorSD(1, dt));
        CHECK(ABS(kmc.getTotProb() - 0.0007299123928392682) < SMALL);
    }
    for (int i = 0; i < 4; ++i) {
        delete ep_j[i];
    }
}

/**********************************************************************
 *                        Test case scenarios                         *
 **********************************************************************/

TEST_CASE("3 crossing perpendicular rods with protein in center",
          "[three_perp_rods]") {
    // Arrange
    double dt = .0001;
    double KBT = 1.;
    ExampleRod *ep_j[3];
    for (int i = 0; i < 3; ++i) {
        ep_j[i] = new ExampleRod;
        *(ep_j[i]) = MockRod<ExampleRod>(i);
        ep_j[i]->length = 20;
    }
    ep_j[0]->direction[0] = 1;
    ep_j[0]->direction[1] = 0;
    ep_j[0]->direction[2] = 0;

    ep_j[1]->direction[0] = 0;
    ep_j[1]->direction[1] = 1;
    ep_j[1]->direction[2] = 0;

    ep_j[2]->direction[0] = 0;
    ep_j[2]->direction[1] = 0;
    ep_j[2]->direction[2] = 1;

    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * ep_j[0]->radius);
    xlink.LUTablePtr = &LUT;
    std::vector<int> Uniquefilter{1, 1, 1};
    // Act
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 3, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        for (int i = 0; i < 3; ++i) {
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(i) - 0.0001909859317102744) < SMALL);
            // Check minimum distance
            CHECK(ABS(kmc.getDistMin(i) - 0.0) < SMALL);
            // Check mu distances along rods
            CHECK(ABS(kmc.getMu(i) - 0.0) < SMALL);
        }
        // Check total probability
        REQUIRE(ABS(kmc.getTotProb() - 0.0005729577951308232) < SMALL);
        // Check binding to rod 0
        double roll = .00001;
        double bindpos;
        int rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 0);
        CHECK(ABS(bindpos + 0.4476401224401701) < SMALL);
        // Check binding to rod 1
        roll = .0003;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 1);
        CHECK(ABS(bindpos - 0.07079632679489645) < SMALL);
        // Check out of range binding returns -1
        roll = .9;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == -1);
    }
    SECTION("1 head bound to first rod") {
        int end_bound = 0;
        double end_loc = 0;
        // Bind to center of first rod
        xlink.setBind(end_bound, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos,
                      0, ep_j[0]->length, ep_j[0]->rank);
        SECTION("1->0 unbinding") {
            double rollVec[] = {.1, .125, .125};
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 3,
                                xlink.getRcutUS());
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorDS(end_bound, dt, KBT));
            CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
            // Ignore rescaling since this is a test
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos[] = {0.1085452196603419, 0.1085452196603419,
                                  -0.1740595812604792};
            CHECK(ABS(pos[0] - check_pos[0]) < SMALL);
            CHECK(ABS(pos[1] - check_pos[1]) < SMALL);
            CHECK(ABS(pos[2] - check_pos[2]) < SMALL);
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos2[] = {0., 0., -0.232079441680639};
            CHECK(ABS(pos[0] - check_pos2[0]) < SMALL);
            CHECK(ABS(pos[1] - check_pos2[1]) < SMALL);
            CHECK(ABS(pos[2] - check_pos2[2]) < SMALL);
        }
        SECTION("1->2 binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 3,
                                xlink.getRcutSD());
            kmc.CalcTotProbsSD(ep_j, Uniquefilter, xlink.idBind[end_bound],
                               xlink.lambda, xlink.kappa, 1. / KBT,
                               xlink.freeLength,
                               xlink.getBindingFactorSD(end_bound, dt));
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(0)) == 0.0); // No binding, already attached
            CHECK(ABS(kmc.getProbs(1) - 0.0004899204301239162) < SMALL);
            CHECK(ABS(kmc.getProbs(2) - 0.0004899204301239162) < SMALL);
            // Check total probability
            CHECK(ABS(kmc.getTotProb() - 0.0009798408602478324) < SMALL);
            // Check minimum distances
            CHECK(ABS(kmc.getDistMin(1) - 0.0) < SMALL);
            CHECK(ABS(kmc.getDistMin(2) - 0.0) < SMALL);
            // Check mu distances along rods
            CHECK(ABS(kmc.getMu(1) - 0.0) < SMALL);
            CHECK(ABS(kmc.getMu(2) - 0.0) < SMALL);
            // Choose a (random) vector
            double rollVec[] = {.0001, .125, .125};
            // TODO: Implement LUT or taylor expansion <31-12-18, ARL>
            // Check placement of head Check ID of rod
        }
    }
    SECTION("2 heads bound to first and second rod") {
        // Bind to center of first rod
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        double rollVec[] = {.0001, .125, .125};
        xlink.setBind(1, ep_j[1]->gid, ep_j[1]->direction, ep_j[1]->pos, 0,
                      ep_j[1]->length, ep_j[1]->rank);
        // Check individual probabilities
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 3, xlink.getRcutSD());
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, dt, KBT));
        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(ABS(err[0]) < SMALL);
        CHECK(ABS(err[1]) < SMALL);
        CHECK(ABS(err[2]) < SMALL);
    }
    for (int i = 0; i < 3; ++i) {
        delete ep_j[i];
    }
}

TEST_CASE("4 parallel rods separated by a rod diameter on the sides with "
          "protein in center",
          "[four_parallel_rods]") {
    // Arrange
    double dt = .0001;
    double KBT = 1.;
    ExampleRod *ep_j[4];
    for (int i = 0; i < 4; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
        ep_j[i]->length = 20;
    }
    // Set position of rods
    ep_j[0]->pos[1] = .25;
    ep_j[1]->pos[1] = -.25;
    ep_j[2]->pos[2] = .25;
    ep_j[3]->pos[2] = -.25;
    // Set directions of rods
    ep_j[0]->direction[0] = 1;
    ep_j[1]->direction[0] = -1; // One anti-parallel
    ep_j[2]->direction[0] = 1;
    ep_j[3]->direction[0] = 1;
    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * ep_j[0]->radius);
    xlink.LUTablePtr = &LUT;
    std::vector<int> Uniquefilter{1, 1, 1, 1};
    // Act
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 4, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        for (int i = 0; i < 4; ++i) {
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(i) - 0.0001653986686265376) < SMALL);
            // Check minimum distance
            CHECK(ABS(kmc.getDistMin(i) - 0.25) < SMALL);
            // Check mu distances along rods
            CHECK(ABS(kmc.getMu(i) - 0.0) < SMALL);
        }
        // Check total probability
        CHECK(ABS(kmc.getTotProb() - (4 * 0.0001653986686265376)) < SMALL);
        // Check binding to rod 0
        double roll = .00005;
        double bindpos;
        int rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 0);
        CHECK(ABS(bindpos + 0.1712133140930698) < SMALL);
        // Check binding to rod 1
        roll = .0002;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 1);
        CHECK(ABS(bindpos + 0.2518405544800601) < SMALL);
        // Check binding to rod 2
        roll = .00035;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 2);
        // Check out of range binding returns -1
        CHECK(ABS(bindpos + 0.3324677948670505) < SMALL);
        roll = .01;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == -1);
    }
    SECTION("1 head bound to first rod") {
        int end_bound = 0;
        double end_loc = 0;
        // Bind to center of first rod
        xlink.setBind(end_bound, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos,
                      0, ep_j[0]->length, ep_j[0]->rank);
        SECTION("1->0 unbinding") {
            double rollVec[] = {.1, .125, .125};
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 4,
                                xlink.getRcutUS());
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorDS(end_bound, dt, KBT));
            CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
            // Check unbinding
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            double *endPos = xlink.posEndBind[end_bound];
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos[] = {0.1085452196603419, 0.3585452196603419,
                                  -0.1740595812604792};
            double err[3] = {};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos[i];
            }
            CHECK(ABS(err[0]) < SMALL);
            CHECK(ABS(err[1]) < SMALL);
            CHECK(ABS(err[2]) < SMALL);
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos1[] = {0., 0.25, -0.232079441680639};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos1[i];
            }
            CHECK(ABS(err[0]) < SMALL);
            CHECK(ABS(err[1]) < SMALL);
            CHECK(ABS(err[2]) < SMALL);
        }
        SECTION("1->2 binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 4,
                                xlink.getRcutSD());
            kmc.CalcTotProbsSD(ep_j, Uniquefilter, xlink.idBind[end_bound],
                               xlink.lambda, xlink.kappa, 1. / KBT,
                               xlink.freeLength,
                               xlink.getBindingFactorSD(end_bound, dt));
            // Check minimum distances
            CHECK(ABS(kmc.getDistMin(1) - .5) < SMALL);
            CHECK(ABS(kmc.getDistMin(2) - 0.3535533905932738) < SMALL);
            CHECK(ABS(kmc.getDistMin(3) - 0.3535533905932738) < SMALL);
            // Check mu distances along rods
            CHECK(ABS(kmc.getMu(1) - 0.0) < SMALL);
            CHECK(ABS(kmc.getMu(2) - 0.0) < SMALL);
            CHECK(ABS(kmc.getMu(3) - 0.0) < SMALL);
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(0)) == 0.0);
            CHECK(ABS(kmc.getProbs(1) - 0.000511748258908916) < SMALL);
            CHECK(ABS(kmc.getProbs(2) - 0.0005021912513360592) < SMALL);
            CHECK(ABS(kmc.getProbs(3) - 0.0005021912513360592) < SMALL);
            // Check total probability
            CHECK(ABS(kmc.getTotProb() - 0.001516130761581034) < SMALL);
            // Choose a (random) vector
            double rollVec[] = {.0001, .125, .125};
        }
    }
    SECTION("2 heads bound to first and second rod") {
        // Bind to center of first rod
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        double rollVec[] = {.0001, .125, .125};
        xlink.setBind(1, ep_j[1]->gid, ep_j[1]->direction, ep_j[1]->pos, 0,
                      ep_j[1]->length, ep_j[1]->rank);
        // Check individual probabilities
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 4, xlink.getRcutSD());
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, dt, KBT));
        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(ABS(err[0]) < SMALL);
        CHECK(ABS(err[1]) < SMALL);
        CHECK(ABS(err[2]) < SMALL);
    }
    for (int i = 0; i < 4; ++i) {
        delete ep_j[i];
    }
}

TEST_CASE("6 perpendicular rods surrounding a sphere of a rod diameter with "
          "protein in center",
          "[six_pointing_perp_rods]") {
    // Arrange
    double dt = .0001;
    double KBT = 1.;
    ExampleRod *ep_j[6];
    for (int i = 0; i < 6; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
        ep_j[i]->length = 20;
    }
    // Set position of rods
    ep_j[0]->pos[0] = 10.25;
    ep_j[0]->pos[1] = 0.0;
    ep_j[0]->pos[2] = 0.0;

    ep_j[1]->pos[0] = -10.25;
    ep_j[1]->pos[1] = 0.0;
    ep_j[1]->pos[2] = 0.0;

    ep_j[2]->pos[0] = 0.0;
    ep_j[2]->pos[1] = 10.25;
    ep_j[2]->pos[2] = 0.0;

    ep_j[3]->pos[0] = 0.0;
    ep_j[3]->pos[1] = -10.25;
    ep_j[3]->pos[2] = 0.0;

    ep_j[4]->pos[0] = 0.0;
    ep_j[4]->pos[1] = 0.0;
    ep_j[4]->pos[2] = 10.25;

    ep_j[5]->pos[0] = 0.0;
    ep_j[5]->pos[1] = 0.0;
    ep_j[5]->pos[2] = -10.25;
    // Set directions of rods
    ep_j[0]->direction[0] = 1.0;
    ep_j[0]->direction[1] = 0.0;
    ep_j[0]->direction[2] = 0.0;

    ep_j[1]->direction[0] = -1.0;
    ep_j[1]->direction[1] = 0.0;
    ep_j[1]->direction[2] = 0.0;

    ep_j[2]->direction[0] = 0.0;
    ep_j[2]->direction[1] = 1.0;
    ep_j[2]->direction[2] = 0.0;

    ep_j[3]->direction[0] = 0.0;
    ep_j[3]->direction[1] = -1.0;
    ep_j[3]->direction[2] = 0.0;

    ep_j[4]->direction[0] = 0.0;
    ep_j[4]->direction[1] = 0.0;
    ep_j[4]->direction[2] = 1.0;

    ep_j[5]->direction[0] = 0.0;
    ep_j[5]->direction[1] = 0.0;
    ep_j[5]->direction[2] = 1.0;
    std::vector<int> Uniquefilter{1, 1, 1, 1, 1, 1};
    // Set protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // ProteinType *pType = &xlink.property;
    LookupTable LUT;
    LUT.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT, xlink.freeLength,
             2. * ep_j[0]->radius);
    xlink.LUTablePtr = &LUT;
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 6, xlink.getRcutUS());
        kmc.CalcTotProbsUS(ep_j, Uniquefilter, xlink.getBindingFactorUS(0, dt));
        for (int i = 0; i < 6; ++i) {
            // Check minimum distance
            CHECK(ABS(kmc.getDistMin(i) - 0.25) < SMALL);
            CHECK(ABS(kmc.getDistPerp(i) - 0.0) < SMALL);
            // Check mu distances along rods
            if (i == 5) {
                CHECK(ABS(kmc.getMu(i) - 10.25) < SMALL);
            } else {
                CHECK(ABS(kmc.getMu(i) + 10.25) < SMALL);
            }
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(i) - 0.0000477464829275686) < SMALL);
        }
        // Check total probability
        CHECK(ABS(kmc.getTotProb() - (6 * 0.0000477464829275686)) < SMALL);
        // Check binding to rod 0
        double roll = .00001;
        double bindpos;
        int rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 0);
        CHECK(ABS((bindpos) + 9.94764012244017) < SMALL);
        // Check binding to rod 1
        roll = .00005;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 1);
        CHECK(ABS(bindpos + 9.98820061220085) < SMALL);
        // Check binding to rod 5
        roll = .00028;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == 5);
        CHECK(ABS(bindpos - 9.966076571675236) < SMALL);
        // Check out of range binding returns -1
        roll = .01;
        rod_i = kmc.whichRodBindUS(ep_j, bindpos, roll);
        CHECK(rod_i == -1);
    }
    SECTION("1 head bound to first rod") {
        int end_bound = 0;
        // Bind to minus end of the first rod
        xlink.setBind(end_bound, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos,
                      ep_j[0]->length * -.5, ep_j[0]->length, ep_j[0]->rank);
        SECTION("1->0 unbinding") {
            double rollVec[] = {.1, .125, .125};
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 6,
                                xlink.getRcutUS());
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorDS(end_bound, dt, KBT));
            CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
            // Check unbinding
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            // const double *endPos = xlink.posEndBind[end_bound];
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos[] = {0.3585452196603419, 0.1085452196603419,
                                  -0.1740595812604792};
            double err[3] = {};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos[i];
            }
            CHECK(ABS(err[0]) < SMALL);
            CHECK(ABS(err[1]) < SMALL);
            CHECK(ABS(err[2]) < SMALL);
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), rollVec, pos);
            double check_pos1[] = {.25, 0., -0.232079441680639};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos1[i];
            }
            CHECK(ABS(err[0]) < SMALL);
            CHECK(ABS(err[1]) < SMALL);
            CHECK(ABS(err[2]) < SMALL);
        }
        SECTION("1->2 binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 6,
                                xlink.getRcutSD());
            // Calculate probabilities
            kmc.CalcTotProbsSD(ep_j, Uniquefilter, xlink.idBind[end_bound],
                               xlink.lambda, xlink.kappa, 1. / KBT,
                               xlink.freeLength,
                               xlink.getBindingFactorSD(end_bound, dt));
            // Check minimum distances
            CHECK(ABS(kmc.getDistMin(1) - .5) < SMALL);
            CHECK(ABS(kmc.getDistMin(2) - 0.3535533905932738) < SMALL);
            CHECK(ABS(kmc.getDistMin(3) - 0.3535533905932738) < SMALL);
            CHECK(ABS(kmc.getDistMin(4) - 0.3535533905932738) < SMALL);
            CHECK(ABS(kmc.getDistMin(5) - 0.3535533905932738) < SMALL);
            // Check mu distances along rods
            CHECK(ABS(kmc.getMu(1) - (ep_j[1]->length * -.5 - .5)) < SMALL);
            CHECK(ABS(kmc.getMu(2) - (ep_j[1]->length * -.5 - .25)) < SMALL);
            CHECK(ABS(kmc.getMu(3) - (ep_j[1]->length * -.5 - .25)) < SMALL);
            CHECK(ABS(kmc.getMu(4) - (ep_j[1]->length * -.5 - .25)) < SMALL);
            CHECK(ABS(kmc.getMu(5) - (ep_j[1]->length * .5 + .25)) < SMALL);
            // Check perpendicular distance
            CHECK(ABS(kmc.getDistPerp(1) - 0.0) < SMALL);
            CHECK(ABS(kmc.getDistPerp(2) - .25) < SMALL);
            CHECK(ABS(kmc.getDistPerp(3) - .25) < SMALL);
            CHECK(ABS(kmc.getDistPerp(4) - .25) < SMALL);
            CHECK(ABS(kmc.getDistPerp(5) - .25) < SMALL);
            // Check individual probabilities
            CHECK(ABS(kmc.getProbs(0)) == 0.0);
            CHECK(ABS(kmc.getProbs(1) - 0.000233916745498158) < SMALL);
            CHECK(ABS(kmc.getProbs(2) - 0.0002425720894371809) < SMALL);
            CHECK(ABS(kmc.getProbs(3) - 0.0002425720894371809) < SMALL);
            CHECK(ABS(kmc.getProbs(4) - 0.0002425720894371809) < SMALL);
            CHECK(ABS(kmc.getProbs(5) - 0.0002425720894371809) < SMALL);
            //// Check total probability
            CHECK(ABS(kmc.getTotProb() - 0.001204205103246882) < SMALL);
            // Choose a (random) vector
            double rollVec[] = {.0001, .125, .125};
        }
    }
    SECTION("2 heads bound to first and second rod") {
        // Bind to minus of first rod
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos,
                      ep_j[0]->length * -.5, ep_j[0]->length, ep_j[0]->rank);
        double rollVec[] = {.0001, .125, .125};
        xlink.setBind(1, ep_j[5]->gid, ep_j[5]->direction, ep_j[5]->pos,
                      ep_j[5]->length * .5, ep_j[5]->length, ep_j[5]->rank);
        // Check individual probabilities
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 6, xlink.getRcutSD());
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, dt, KBT));
        CHECK(ABS(kmc.getTotProb() - .0001) < SMALL);
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(ABS(err[0]) < SMALL);
        CHECK(ABS(err[1]) < SMALL);
        CHECK(ABS(err[2]) < SMALL);
    }
    for (int i = 0; i < 6; ++i) {
        delete ep_j[i];
    }
}

