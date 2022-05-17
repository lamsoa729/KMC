/**********************************************************************
 *                     Unit testing for KMC<ExampleRod> code *
 **********************************************************************/

#include "KMC/ExampleRod.hpp"
#include "KMC/ExampleSphere.hpp"
#include "KMC/helpers.hpp"
#include "KMC/integrals.hpp"
#include "KMC/kmc.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/lut_filler_edep.hpp"
#include "KMC/macros.hpp"
#include "example_objs/ExampleXlink.hpp"

#include "catch.hpp"

#include <array>
#include <cassert>
#include <cmath>

double test_bind_prob(double lm, double sbound0, double sbound1, double lambda,
                      double kappa, double beta, double ell0, double bindFactor,
                      double dt) {
    double prob = integral(lm, sbound0, sbound1,
                           (1. - lambda) * kappa * .5 * beta, ell0) *
                  bindFactor * dt;
    return prob;
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
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);

        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0);
        REQUIRE(prob == Approx(0.0001909859317102744).epsilon(1e-8));
    }
    SECTION(
        "Test binding probability when crosslinker is above and outside rod ") {
        xlink.pos[2] = 2.0;
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);

        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
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
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == rod.length * .5);
        CHECK(kmc.getDistMin(0) == 0);
        CHECK(prob == Approx(0.0000954929658551372).epsilon(1e-8));
    }
    SECTION(
        "Test binding probability when crosslinker is at minus end of rod.") {
        double pos[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
            pos[i] = rod.pos[i] - (rod.direction[i] * rod.length * .5);
        }
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == rod.length * -.5);
        CHECK(kmc.getDistMin(0) == 0);
        REQUIRE(prob == Approx(0.0000954929658551372).epsilon(1e-8));
    }
    SECTION("Test binding probability when crosslinker is on surface of rod "
            "above its center.") {
        double pos[3] = {0, 0, 0.5};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.5);
        REQUIRE(prob == 0);
    }
    SECTION("Test binding probability when crosslinker is on surface of rod's "
            "tip.") {
        double pos[3] = {rod.length * 0.5, 0, 0.5};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == rod.length * .5);
        CHECK(kmc.getDistMin(0) == 0.5);
        REQUIRE(prob == 0);
    }
    SECTION("Test binding probability when crosslinker is between surface of "
            "rod and rod axis above its center.") {
        double pos[3] = {0, 0, 0.25};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.25);
        REQUIRE(prob == Approx(0.0001653986686265376).epsilon(1e-8));
    }
    SECTION("Test binding probability when crosslinker binds to a rod smaller "
            "than cutoff length") {
        double pos[3] = {0, 0, 0};
        rod.length = .5;
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        double prob = kmc.CalcProbRodUS(0, rod, xlink.getBindingFactorUS(0));
        CHECK(kmc.getMu(0) == 0);
        CHECK(kmc.getDistMin(0) == 0.0);
        REQUIRE(prob == Approx(0.0000954929658551372).epsilon(SMALL));
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
        KMC<ExampleRod> kmc(xlink.posEndBind[0], dt);
        kmc.CalcProbSU(xlink.getUnbindingFactorDS(0, KBT));

        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
    }
}

TEST_CASE("Test CalcProbSD for KMC<ExampleRod> class", "[calc_prob_sd]") {
    double dt = .0001;
    double KBT = 1.;
    // Protein data
    ExampleXlink xlink;
    xlink.setMockXlink();
    // Sylinder data
    ExampleRod rod0, rod1;
    rod0 = MockRod<ExampleRod>(0);
    rod0.length = 20;
    rod1 = MockRod<ExampleRod>(1);
    rod1.length = 20;
    // Initialize lookup table
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * rod0.radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;

    SECTION("Test binding to 2 parallel overlapping rods") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        // Check to make sure that probability matches up
        double test_prob = test_bind_prob(
            0, -.5 * rod0.length, .5 * rod0.length, xlink.lambda, xlink.kappa,
            1. / KBT, xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // CHECK(prob == Approx(0.0004899204301239162).epsilon(SMALL));
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 parallel vertically separated rods.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Change the position of rod1 vertically
        double sep = 1.0;
        rod1.pos[0] = 0;
        rod1.pos[1] = sep;
        rod1.pos[2] = 0;
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        // Check to make sure that probability matches up
        double test_prob = test_bind_prob(
            sep, -.5 * rod0.length, .5 * rod0.length, xlink.lambda, xlink.kappa,
            1. / KBT, xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 parallel rods with crosslinker bound at rod "
            "plus end.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Change the position of rod1 vertically
        double sep = 1.0;
        rod1.pos[0] = 0;
        rod1.pos[1] = sep;
        rod1.pos[2] = 0;
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 1, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        double test_prob = test_bind_prob(
            sep, 0, .5 * rod0.length, xlink.lambda, xlink.kappa, 1. / KBT,
            xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // Check to make sure that probability matches up
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 parallel rods with crosslinker bound at rod "
            "minus end.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = -rod0.length * .5;
        // Change the position of rod1 vertically
        double sep = 1.0;
        rod1.pos[0] = 0;
        rod1.pos[1] = sep;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        std::cout << xlink.getRcutSD() << std::endl;
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        double test_prob = test_bind_prob(
            sep, 0, .5 * rod0.length, xlink.lambda, xlink.kappa, 1. / KBT,
            xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // Check to make sure that probability matches up
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 parallel inline rods with crosslinker bound at "
            "rod plus end and other rod shifted by rod length + D.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        // Change the position of rod1 horizontally
        double sep = 1.0;
        rod1.pos[0] = rod0.length + sep;
        rod1.pos[1] = 0;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        double test_prob = test_bind_prob(
            0, sep, .5 * rod0.length, xlink.lambda, xlink.kappa, 1. / KBT,
            xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // Check to make sure that probability matches up
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 parallel with crosslinker bound at "
            "rod plus end and other rod shifted horizontially by rod length + "
            "D and vertically by D") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = rod0.length * .5;
        // Change the position of rod1 horizontally
        double sep = 1.0;
        rod1.pos[0] = rod0.length + sep;
        rod1.pos[1] = sep;
        rod1.pos[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        double test_prob = test_bind_prob(
            sep, sep, .5 * rod0.length, xlink.lambda, xlink.kappa, 1. / KBT,
            xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // Check to make sure that probability matches up
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION("Test binding to 2 perpendicular rods crosslinked in the center "
            "and one rod shifted over by D.") {
        // Bind protein head 0 to the center of rod0
        int end_bound = 0;
        double end_loc = 0;
        // Change the position and direction of rod0
        double sep = 1.0;
        rod0.pos[0] = 0;
        rod0.pos[1] = sep;
        rod0.pos[2] = 0;
        rod0.direction[0] = 0;
        rod0.direction[1] = 1.0;
        rod0.direction[2] = 0;
        xlink.setBind(end_bound, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        // Create KMC<ExampleRod> object for protein
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        // Calculate binding of protein head 1 to rod 1
        double prob =
            kmc.CalcProbRodSD(0, rod1, xlink.lambda, xlink.kappa, 1. / KBT,
                              xlink.freeLength, xlink.getBindingFactorSD(1));
        double test_prob = test_bind_prob(
            sep, -.5 * rod0.length, .5 * rod0.length, xlink.lambda, xlink.kappa,
            1. / KBT, xlink.freeLength, xlink.getBindingFactorSD(1), dt);
        // Check to make sure that probability matches up
        CHECK(prob == Approx(test_prob).epsilon(SMALL));
    }
    SECTION(
        "Exception handling when no lookup table is given for SD binding.") {
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 0, xlink.getRcutSD(), dt);
        REQUIRE_THROWS(
            kmc.LUCalcProbRodSD(0, rod0, xlink.getBindingFactorUS(0)));
        REQUIRE_THROWS(kmc.RandomBindPosSD(0, .1));
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
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * rod0.radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;

    SECTION("Test unbinding of end 0 from rod0.") {
        double end_loc = 0;
        xlink.setBind(0, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        xlink.setBind(1, rod1.gid, rod1.direction, rod1.pos, end_loc,
                      rod0.length, rod0.rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(0, KBT));
        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
    }
    SECTION("Test unbinding of end 1 from rod1.") {
        double end_loc = 0;
        xlink.setBind(0, rod0.gid, rod0.direction, rod0.pos, end_loc,
                      rod0.length, rod0.rank);
        xlink.setBind(1, rod1.gid, rod1.direction, rod1.pos, end_loc,
                      rod0.length, rod0.rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 2, xlink.getRcutSD(), dt);
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, KBT));
        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
    }
}

TEST_CASE("Test binding probabilities of one head when near multiple rods.",
          "[calc01tot_probs]") {
    double dt = .0001;
    double KBT = 1.;
    std::vector<ExampleRod *> ep_j(4);
    for (int i = 0; i < 4; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
    }
    // Initialize crosslink object
    ExampleXlink xlink;
    xlink.setMockXlink();
    // Intialize lookup table for binding calculations
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(xlink.getExponentFactor() / KBT, xlink.freeLength,
                    2. * ep_j[0]->radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;
    // Initialize other vectors for calculations
    std::vector<double> bindFactors(4, xlink.getBindingFactorUS(0));
    // Testing
    SECTION("Two semi-overlapping parallel rods with protein in between.") {
        double sep = .25;
        ep_j[0]->pos[0] = 0.0;
        ep_j[0]->pos[1] = 2. * sep;
        ep_j[0]->pos[2] = 0.0;
        double pos[3] = {0, sep, 0};

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 2, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        kmc.CalcTotProbsUS(cep_j, bindFactors);
        REQUIRE(kmc.getTotProb() ==
                Approx(0.000330797337253075).epsilon(SMALL));
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

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        double pos[3] = {0, 0.25, 0};
        xlink.setUnBindPos(pos);
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 4, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        kmc.CalcTotProbsUS(cep_j, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.0006615946745061504).epsilon(SMALL));
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
    // Initialize rods
    std::vector<ExampleRod *> ep_j(4);
    // ExampleRod *ep_j[4];
    for (int i = 0; i < 4; ++i) {
        ep_j[i] = new ExampleRod;
        *ep_j[i] = MockRod<ExampleRod>(i);
    }
    // Initialize crosslinks
    ExampleXlink xlink;
    xlink.setMockXlink();

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * ep_j[0]->radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;

    std::vector<double> bindFactors(4, xlink.getBindingFactorUS(1));

    SECTION("No rods nearby") {
        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 1, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, cep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0).margin(SMALL));
    }
    SECTION("One rod nearby and above bound crosslinker") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 2, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, cep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.0005481166497130916).epsilon(SMALL));
    }
    SECTION("Two rods nearby, above and below bound crosslinker") {
        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, ep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.001096233299426183).epsilon(SMALL));
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

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        // Bind crosslinker to rod
        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, cep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.001096233299426183).epsilon(SMALL));
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
        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, cep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
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

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        KMC<ExampleRod> kmc(xlink.posEndBind[0], 3, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, ep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);

        CHECK(kmc.getTotProb() == Approx(0).margin(SMALL));
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

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);

        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, ep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.001644349949139275).epsilon(SMALL));
    }
    SECTION("Three rods nearby above, below, and on, all inside rcut.") {
        ep_j[1]->length = ep_j[2]->length = ep_j[3]->length = 2.0;

        ep_j[1]->pos[0] = 0;
        ep_j[1]->pos[1] = 1.0;
        ep_j[1]->pos[2] = 0;

        ep_j[2]->pos[0] = 0;
        ep_j[2]->pos[1] = -1.0;
        ep_j[2]->pos[2] = 0;

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        xlink.setBind(0, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos, 0,
                      ep_j[0]->length, ep_j[0]->rank);
        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, ep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.0003456279751132963).epsilon(SMALL));
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

        // Copy ep_j ptr to const ptr vector so it can run with KMC
        std::vector<const ExampleRod *> cep_j;
        cep_j.assign(ep_j.begin(), ep_j.end());

        KMC<ExampleRod> kmc(xlink.posEndBind[0], 4, xlink.getRcutSD(), dt);

        kmc.CalcTotProbsSD(cep_j, ep_j[0]->gid, xlink.lambda, xlink.kappa,
                           1. / KBT, xlink.freeLength, bindFactors);
        CHECK(kmc.getTotProb() == Approx(0.0007299123928392682).epsilon(SMALL));
    }
    // Cleanup
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
    std::vector<ExampleRod *> ep_j(3);
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

    // Copy ep_j ptr to const ptr vector so it can run with KMC
    std::vector<const ExampleRod *> cep_j;
    cep_j.assign(ep_j.begin(), ep_j.end());

    ExampleXlink xlink;
    xlink.setMockXlink();

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * ep_j[0]->radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;

    // Act
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 3, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        std::vector<double> bindFactorsUS(3, xlink.getBindingFactorUS(0));
        kmc.CalcTotProbsUS(cep_j, bindFactorsUS);
        for (int i = 0; i < 3; ++i) {
            // Check individual probabilities
            CHECK(kmc.getProbs(i) ==
                  Approx(0.0001909859317102744).epsilon(SMALL));
            // Check minimum distance
            CHECK(kmc.getDistMin(i) == Approx(0.0).epsilon(SMALL));
            // Check mu distances along rods
            CHECK(kmc.getMu(i) == Approx(0.0).epsilon(SMALL));
        }
        // Check total probability
        double totProb = kmc.getTotProb();
        REQUIRE(totProb == Approx(0.0005729577951308232).epsilon(SMALL));
        // Check binding to rods
        double roll = 0;
        double prevRodProb = 0;
        for (int i = 0; i < 3; ++i) {
            double RodProb = kmc.getProbs(i);
            for (double j = 0; j < 1.; j += .2) {
                double roll = ((j * RodProb) + prevRodProb) / totProb;
                double calcBindPos = 2. * xlink.rc * (j - .5);
                double bindpos;
                int rod_i = kmc.whichObjBindUS(cep_j, bindpos, roll);
                CHECK(rod_i == i);
                CHECK(bindpos == Approx(calcBindPos).epsilon(SMALL));
            }
            prevRodProb += RodProb;
        }
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
                                xlink.getRcutUS(), dt);
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorDS(end_bound, KBT));
            CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
            // Ignore rescaling since this is a test
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos[] = {0.1085452196603419, 0.1085452196603419,
                                  -0.1740595812604792};
            CHECK(pos[0] == Approx(check_pos[0]).epsilon(SMALL));
            CHECK(pos[1] == Approx(check_pos[1]).epsilon(SMALL));
            CHECK(pos[2] == Approx(check_pos[2]).epsilon(SMALL));
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos2[] = {0., 0., -0.232079441680639};
            CHECK(pos[0] == Approx(check_pos2[0]).margin(SMALL));
            CHECK(pos[1] == Approx(check_pos2[1]).margin(SMALL));
            CHECK(pos[2] == Approx(check_pos2[2]).margin(SMALL));
        }
        SECTION("1->2 binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 3,
                                xlink.getRcutSD(), dt);
            std::vector<double> bindFactorsSD(
                3, xlink.getBindingFactorSD(1 - end_bound));
            kmc.CalcTotProbsSD(cep_j, xlink.idBind[end_bound], xlink.lambda,
                               xlink.kappa, 1. / KBT, xlink.freeLength,
                               bindFactorsSD);
            // Check individual probabilities
            CHECK(kmc.getProbs(0) == 0.0); // No binding, already attached
            CHECK(kmc.getProbs(1) ==
                  Approx(0.0004899204301239162).epsilon(SMALL));
            CHECK(kmc.getProbs(2) ==
                  Approx(0.0004899204301239162).epsilon(SMALL));
            // Check total probability
            CHECK(kmc.getTotProb() ==
                  Approx(0.0009798408602478324).epsilon(SMALL));
            // Check minimum distances
            CHECK(kmc.getDistMin(1) == Approx(0.0).epsilon(SMALL));
            CHECK(kmc.getDistMin(2) == Approx(0.0).epsilon(SMALL));
            // Check mu distances along rods
            CHECK(kmc.getMu(1) == Approx(0.0).epsilon(SMALL));
            CHECK(kmc.getMu(2) == Approx(0.0).epsilon(SMALL));
            // Choose a (random) vector
            double rollVec[] = {.0001, .125, .125};
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
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 3, xlink.getRcutSD(), dt);
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, KBT));
        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(err[0] == Approx(0).margin(SMALL));
        CHECK(err[1] == Approx(0).margin(SMALL));
        CHECK(err[2] == Approx(0).margin(SMALL));
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
    std::vector<ExampleRod *> ep_j(4);
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

    // Copy ep_j ptr to const ptr vector so it can run with KMC
    std::vector<const ExampleRod *> cep_j;
    cep_j.assign(ep_j.begin(), ep_j.end());

    ExampleXlink xlink;
    xlink.setMockXlink();
    // Initialize lookup table
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * ep_j[0]->radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;

    // Act
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 4, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        std::vector<double> bindFactorsUS(4, xlink.getBindingFactorUS(0));
        kmc.CalcTotProbsUS(cep_j, bindFactorsUS);
        for (int i = 0; i < 4; ++i) {
            // Check individual probabilities
            REQUIRE(kmc.getProbs(i) ==
                    Approx(0.0001653986686265376).epsilon(SMALL));
            // Check minimum distance
            CHECK(kmc.getDistMin(i) == Approx(0.25).epsilon(SMALL));
            // Check mu distances along rods
            CHECK(kmc.getMu(i) == Approx(0.0).epsilon(SMALL));
        }
        // Check total probability
        double totProb = kmc.getTotProb();
        REQUIRE(totProb == Approx(4. * kmc.getProbs(0)).epsilon(SMALL));
        // Check binding to rods
        double prevRodProb = 0;
        for (int i = 0; i < 3; ++i) {
            double modRC = sqrt(SQR(xlink.rc) - SQR(kmc.getDistMin(i)));
            double RodProb = kmc.getProbs(i);
            for (double j = 0; j < 1.; j += .2) {
                double roll = ((j * RodProb) + prevRodProb) / totProb;
                double calcBindPos = 2. * modRC * (j - .5);
                double bindpos;
                int rod_i = kmc.whichObjBindUS(cep_j, bindpos, roll);
                CHECK(rod_i == i);
                CHECK(bindpos == Approx(calcBindPos).epsilon(SMALL));
            }
            prevRodProb += RodProb;
        }
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
                                xlink.getRcutUS(), dt);
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorSU(end_bound));
            CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
            // Check unbinding
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            double *endPos = xlink.posEndBind[end_bound];
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos[] = {0.1085452196603419, 0.3585452196603419,
                                  -0.1740595812604792};
            double err[3] = {};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos[i];
            }
            CHECK(err[0] == Approx(0).margin(SMALL));
            CHECK(err[1] == Approx(0).margin(SMALL));
            CHECK(err[2] == Approx(0).margin(SMALL));
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos1[] = {0., 0.25, -0.232079441680639};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos1[i];
            }
            CHECK(err[0] == Approx(0).margin(SMALL));
            CHECK(err[1] == Approx(0).margin(SMALL));
            CHECK(err[2] == Approx(0).margin(SMALL));
        }
        SECTION("1->2 binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 4,
                                xlink.getRcutSD(), dt);
            std::vector<double> bindFactorsSD(
                4, xlink.getBindingFactorSD(1 - end_bound));
            kmc.CalcTotProbsSD(cep_j, xlink.idBind[end_bound], xlink.lambda,
                               xlink.kappa, 1. / KBT, xlink.freeLength,
                               bindFactorsSD);
            // Check minimum distances
            CHECK(kmc.getDistMin(1) == Approx(.5).epsilon(SMALL));
            CHECK(kmc.getDistMin(2) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            CHECK(kmc.getDistMin(3) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            // Check mu distances along rods
            CHECK(kmc.getMu(1) == Approx(0.0).margin(SMALL));
            CHECK(kmc.getMu(2) == Approx(0.0).margin(SMALL));
            CHECK(kmc.getMu(3) == Approx(0.0).margin(SMALL));
            // Check individual probabilities
            CHECK(kmc.getProbs(0) == 0.0);
            CHECK(kmc.getProbs(1) ==
                  Approx(0.000511748258908916).epsilon(SMALL));
            CHECK(kmc.getProbs(2) ==
                  Approx(0.0005021912513360592).epsilon(SMALL));
            CHECK(kmc.getProbs(3) ==
                  Approx(0.0005021912513360592).epsilon(SMALL));
            // Check total probability
            CHECK(kmc.getTotProb() ==
                  Approx(0.001516130761581034).epsilon(SMALL));
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
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 4, xlink.getRcutSD(), dt);
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, KBT));
        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(err[0] == Approx(0).margin(SMALL));
        CHECK(err[1] == Approx(0).margin(SMALL));
        CHECK(err[2] == Approx(0).margin(SMALL));
    }
    for (int i = 0; i < 4; ++i) {
        delete ep_j[i];
    }
}

TEST_CASE("6 perpendicular rods surrounding a sphere of a rod radius with "
          "protein in center",
          "[six_pointing_perp_rods]") {
    // Arrange
    double dt = .0001;
    double KBT = 1.;
    // ExampleRod *ep_j[6];
    std::vector<ExampleRod *> ep_j(6);
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

    // Copy ep_j ptr to const ptr vector so it can run with KMC
    std::vector<const ExampleRod *> cep_j;
    cep_j.assign(ep_j.begin(), ep_j.end());

    // Set protein data
    ExampleXlink xlink;
    xlink.setMockXlink();

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init((1. - xlink.lambda) * .5 * xlink.kappa / KBT,
                    xlink.freeLength, 2. * ep_j[0]->radius);
    LookupTable LUT(&lut_filler);
    xlink.LUTablePtr = &LUT;
    SECTION("Unbound") {
        KMC<ExampleRod> kmc(xlink.getPosPtr(), 6, xlink.getRcutUS(),
                            xlink.getDiffU(), dt);
        std::vector<double> bindFactorsUS(6, xlink.getBindingFactorUS(0));
        kmc.CalcTotProbsUS(cep_j, bindFactorsUS);
        for (int i = 0; i < 6; ++i) {
            // Check minimum distance
            REQUIRE(kmc.getDistMin(i) == Approx(0.25).epsilon(SMALL));
            REQUIRE(kmc.getDistPerp(i) == Approx(0.0).epsilon(SMALL));
            // Check mu distances along rods
            if (i == 5) {
                REQUIRE(kmc.getMu(i) == Approx(10.25).epsilon(SMALL));
            } else {
                REQUIRE(kmc.getMu(i) == Approx(-10.25).epsilon(SMALL));
            }
            // Check individual probabilities
            CHECK(kmc.getProbs(i) ==
                  Approx(0.0000477464829275686).epsilon(SMALL));
        }
        double totProb = kmc.getTotProb();
        REQUIRE(totProb == Approx(6. * kmc.getProbs(0)).epsilon(SMALL));
        // Check binding to rods
        double prevRodProb = 0;
        for (int i = 0; i < 6; ++i) {
            double modRC = xlink.rc - kmc.getDistMin(i);
            double RodProb = kmc.getProbs(i);
            for (double j = 0; j < 1.; j += .2) {
                double roll = ((j * RodProb) + prevRodProb) / totProb;
                double calcBindPos = 0;
                if (i == 5) {
                    calcBindPos =
                        kmc.getMu(i) - kmc.getDistMin(i) - (modRC * (1. - j));
                } else {
                    calcBindPos =
                        kmc.getMu(i) + kmc.getDistMin(i) + (modRC * j);
                }
                double bindpos;
                int rod_i = kmc.whichObjBindUS(cep_j, bindpos, roll);
                CHECK(rod_i == i);
                CHECK(bindpos == Approx(calcBindPos).epsilon(SMALL));
            }
            prevRodProb += RodProb;
        }
    }
    SECTION("1 head bound to first rod") {
        int end_bound = 0;
        // Bind to minus end of the first rod
        xlink.setBind(end_bound, ep_j[0]->gid, ep_j[0]->direction, ep_j[0]->pos,
                      ep_j[0]->length * -.5, ep_j[0]->length, ep_j[0]->rank);
        SECTION("S->U unbinding") {
            double rollVec[] = {.1, .125, .125};
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 6,
                                xlink.getRcutUS(), dt);
            // Check individual probabilities
            kmc.CalcProbSU(xlink.getUnbindingFactorDS(end_bound, KBT));
            CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
            // Check unbinding
            // rollVec[0] = (kmc.getTotProb() - rollVec[0]) /
            kmc.getTotProb();
            // Check unbound location
            double pos[3] = {};
            // const double *endPos = xlink.posEndBind[end_bound];
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos[] = {0.3585452196603419, 0.1085452196603419,
                                  -0.1740595812604792};
            double err[3] = {};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos[i];
            }
            CHECK(err[0] == Approx(0).margin(SMALL));
            CHECK(err[1] == Approx(0).margin(SMALL));
            CHECK(err[2] == Approx(0).margin(SMALL));
            // Test 0 edge cases
            rollVec[1] = rollVec[2] = 0;
            kmc.whereUnbindSU(xlink.getRcutUS(), xlink.getDiffU(), rollVec,
                              pos);
            double check_pos1[] = {.25, 0., -0.232079441680639};
            for (int i = 0; i < 3; ++i) {
                err[i] = pos[i] - check_pos1[i];
            }
            CHECK(err[0] == Approx(0).margin(SMALL));
            CHECK(err[1] == Approx(0).margin(SMALL));
            CHECK(err[2] == Approx(0).margin(SMALL));
        }
        SECTION("S->D binding") {
            KMC<ExampleRod> kmc(xlink.posEndBind[end_bound], 6,
                                xlink.getRcutSD(), dt);
            std::vector<double> bindFactorsSD(
                6, xlink.getBindingFactorSD(1 - end_bound));
            // Calculate probabilities
            kmc.CalcTotProbsSD(cep_j, xlink.idBind[end_bound], xlink.lambda,
                               xlink.kappa, 1. / KBT, xlink.freeLength,
                               bindFactorsSD);
            // Check minimum distances
            CHECK(kmc.getDistMin(1) == Approx(.5).epsilon(SMALL));
            CHECK(kmc.getDistMin(2) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            CHECK(kmc.getDistMin(3) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            CHECK(kmc.getDistMin(4) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            CHECK(kmc.getDistMin(5) ==
                  Approx(0.3535533905932738).epsilon(SMALL));
            // Check mu distances along rods
            CHECK(kmc.getMu(1) ==
                  Approx(ep_j[1]->length * -.5 - .5).epsilon(SMALL));
            CHECK(kmc.getMu(2) ==
                  Approx(ep_j[1]->length * -.5 - .25).epsilon(SMALL));
            CHECK(kmc.getMu(3) ==
                  Approx(ep_j[1]->length * -.5 - .25).epsilon(SMALL));
            CHECK(kmc.getMu(4) ==
                  Approx(ep_j[1]->length * -.5 - .25).epsilon(SMALL));
            CHECK(kmc.getMu(5) ==
                  Approx(ep_j[1]->length * .5 + .25).epsilon(SMALL));
            // Check perpendicular distance
            CHECK(kmc.getDistPerp(1) == Approx(0.0).epsilon(SMALL));
            CHECK(kmc.getDistPerp(2) == Approx(.25).epsilon(SMALL));
            CHECK(kmc.getDistPerp(3) == Approx(.25).epsilon(SMALL));
            CHECK(kmc.getDistPerp(4) == Approx(.25).epsilon(SMALL));
            CHECK(kmc.getDistPerp(5) == Approx(.25).epsilon(SMALL));
            // Check individual probabilities
            CHECK(kmc.getProbs(0) == 0.0);
            CHECK(kmc.getProbs(1) ==
                  Approx(0.000233916745498158).epsilon(SMALL));
            CHECK(kmc.getProbs(2) ==
                  Approx(0.0002425720894371809).epsilon(SMALL));
            CHECK(kmc.getProbs(3) ==
                  Approx(0.0002425720894371809).epsilon(SMALL));
            CHECK(kmc.getProbs(4) ==
                  Approx(0.0002425720894371809).epsilon(SMALL));
            CHECK(kmc.getProbs(5) ==
                  Approx(0.0002425720894371809).epsilon(SMALL));
            //// Check total probability
            CHECK(kmc.getTotProb() ==
                  Approx(0.001204205103246882).epsilon(SMALL));
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
        KMC<ExampleRod> kmc(xlink.posEndBind[1], 6, xlink.getRcutSD(), dt);
        kmc.CalcProbDS(xlink.getUnbindingFactorDS(1, KBT));
        CHECK(kmc.getTotProb() == Approx(.0001).epsilon(SMALL));
        double pos[3] = {};
        kmc.whereUnbindDS(xlink.posEndBind[0], pos);
        double err[3] = {};
        for (int i = 0; i < 3; ++i) {
            err[i] = pos[i] - xlink.posEndBind[0][i];
        }
        CHECK(err[0] == Approx(0).margin(SMALL));
        CHECK(err[1] == Approx(0).margin(SMALL));
        CHECK(err[2] == Approx(0).margin(SMALL));
    }
    for (int i = 0; i < 6; ++i) {
        delete ep_j[i];
    }
}

TEST_CASE("Test calculation of U<->S binding radius based on diffusion.",
          "[cutoff_radius]") {
    // Define constants for KMC initialization
    double dt = .001;
    double pos[3] = {0, 0, 0};
    double rc = .5;
    int Npj = 1;
    int n_trials = 10;

    SECTION("Diffusion radius is larger than given capture radius.") {
        for (int i = 0; i < n_trials; ++i) {
            double rad_capture = (i + 1.) * rc;
            double calc_diff = SQR(rad_capture) / (dt * 6.);
            KMC<ExampleRod> kmc_U(pos, 1, rc, calc_diff, dt);
            CHECK(rad_capture ==
                  Approx(kmc_U.getDiffRadius(calc_diff)).epsilon(SMALL));
            CHECK(rad_capture == Approx(kmc_U.getRcutoff()).epsilon(SMALL));
        }
    }
    SECTION("Diffusion radius is smaller than capture radius.") {
        for (int i = 0; i < n_trials; ++i) {
            double diffConst = i;
            // Make radius capture always greater than diffusion radius
            double rad_capture = (6. * diffConst * dt) + 1.;
            KMC<ExampleRod> kmc_U(pos, 1, rad_capture, diffConst, dt);
            REQUIRE(rad_capture > kmc_U.getDiffRadius(diffConst));
            CHECK(rad_capture == kmc_U.getRcutoff());
        }
    }
    SECTION("Edge cases") {
        // No diffusion, should just be rc
        KMC<ExampleRod> kmc0(pos, 1, rc, 0, dt);
        CHECK(rc == kmc0.getRcutoff());

        // Diff radius is 1 but original rc = 0
        double diffRad = 1.;
        double diffConst = 1. / (dt * 6.);
        KMC<ExampleRod> kmc1(pos, 1, 0.0, diffConst, dt);
        CHECK(diffRad == kmc1.getRcutoff());
    }
}

TEST_CASE("Works with pointlike objects.") {
    double rc = .4;
    double dt = .001;
    double pos[3] = {0.1, 0.0, 0.1};
    double bind_vol = M_PI * rc * rc * rc / 0.75;
    double bind_fac = 0.1;
    double lambda = 0.334;
    double kappa = 0.823;
    double beta = 0.610;
    double l0 = 0.08;
    ExampleSphere sphere = MockSphere<ExampleSphere>(0);
    ExampleSphere sphere2 = MockSphere<ExampleSphere>(1);
    std::vector<const ExampleSphere *> spheres = {&sphere, &sphere2};
    std::vector<const ExampleRod *> rods; // empty vector
    std::vector<double> bind_factors = {0.1, 0.1};
    SECTION("KMC constructed with spheres") {
        KMC<ExampleRod, ExampleSphere> kmc(pos, 0, 1, rc, dt);

        // Check to see that it initialized by checking rcutoff
        CHECK(rc == kmc.getRcutoff());
    }
    SECTION("Unbound to Single Probability") {
        KMC<ExampleRod, ExampleSphere> kmc(pos, 0, 1, rc, 1, dt);
        double prob = kmc.CalcProbSphereUS(0, sphere, bind_fac);
        double exp_prob =
            bind_fac * 4.0 * M_PI * SQR(sphere.radius) * dt / bind_vol;

        // Should be rc because diffusion radius < rcutoff
        CHECK(rc == kmc.getRcutoff());
        CHECK(prob == Approx(exp_prob).epsilon(1e-8));
    }
    SECTION("Single to Double Probabilities (direct)") {
        KMC<ExampleRod, ExampleSphere> kmc(pos, 0, 2, rc, 1, dt);
        double expon =
            exp(-0.5 * (1 - lambda) * kappa * beta *
                SQR(sqrt(SQR(pos[0]) + SQR(pos[1]) + SQR(pos[2])) - l0));
        double exp_prob =
            bind_fac * 4.0 * M_PI * SQR(sphere.radius) * expon * dt;
        double prob =
            kmc.CalcProbSphereSD(0, sphere, lambda, kappa, beta, l0, bind_fac);
        kmc.CalcTotProbsSD(rods, spheres, 552, lambda, kappa, beta, l0,
                           bind_factors);

        // Direct calculations
        CHECK(prob == Approx(exp_prob).epsilon(1e-8));
        // just test that total prob works if we have two identical spheres in
        // our vector. could make this a little more comprehensive if needed.
        CHECK(kmc.getTotProb() == Approx(2 * exp_prob).epsilon(1e-8));

        // Lookup table calculations
        LUTFillerAsym lut(256, 256);
        double k_compress = 1.3;
        double f_dep_l = 0.05;

        lut.Init(k_compress, beta * kappa * 0.5, lambda, f_dep_l, l0, 1.0);
        LookupTable LUT(&lut);
        expon = exp(
            -0.5 * kappa * beta *
            ((1 - lambda) *
                 SQR(sqrt(SQR(pos[0]) + SQR(pos[1]) + SQR(pos[2])) - l0) -
             f_dep_l * (sqrt(SQR(pos[0]) + SQR(pos[1]) + SQR(pos[2])) - l0)));
        exp_prob = bind_fac * 4.0 * M_PI * SQR(sphere.radius) * expon * dt;
        KMC<ExampleRod, ExampleSphere> kmc_lut(pos, 0, 2, dt, &LUT);
        prob = kmc_lut.LUCalcProbSphereSD(0, sphere, bind_fac);
        CHECK(prob == Approx(exp_prob).epsilon(1e-8));

        // Again, test total prob for 2 identical spheres
        kmc_lut.LUCalcTotProbsSD(rods, spheres, 300, bind_factors);
        CHECK(kmc_lut.getTotProb() == Approx(2 * exp_prob).epsilon(1e-8));
    }
}

