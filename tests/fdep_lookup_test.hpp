#include "KMC/integrals.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/lut_filler_fdep.hpp"
#include "catch.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include <boost/math/quadrature/gauss_kronrod.hpp>

TEST_CASE("Fdep lookup table filler upper bound test", "[upper_bound]") {
    // Assemble
    const double D = 0.024;
    const double e_fact = .5;
    const double fdep_length = .01;
    const double freelength = 0.05;

    LUTFillerFdep lut_filler(256, 256);
    SECTION("Test SOFT spring") {
        const double M = 0.1 / 0.00411;
        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        REQUIRE(D * lut_filler.getUpperBound() ==
                Approx(1.300682719937549).epsilon(1e-8));
    }
    SECTION("Test MEDIUM spring") {
        const double M = 1. / 0.00411;
        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        REQUIRE(D * lut_filler.getUpperBound() ==
                Approx(0.4596382883076154).epsilon(1e-8));
    }
    SECTION("Test STIFF spring") {
        const double M = 10. / 0.00411;
        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        REQUIRE(D * lut_filler.getUpperBound() ==
                Approx(0.1946667540747286).epsilon(1e-8));
    }
}

TEST_CASE("Fdep lookup table Lookup method test ", "[lookup]") {

    const double D = 0.024;
    const double e_fact = .5;
    const double fdep_length = .01;
    const double freelength = 0.05;
    LUTFillerFdep lut_filler(256, 256);

    double distPerp = 0;
    SECTION("Test SOFT spring") {
        constexpr double errTol = 1e-4;
        const double M = 0.1 / 0.00411;

        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        LookupTable LUT(&lut_filler);
        const double lUB = LUT.getLUCutoff();

        distPerp = 0.04;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
            // Approx(D * fdep_integral(distPerp / D, 0, fact * lUB / D,
            // M * D * D, e_fact, fdep_length / D,
            // freelength / D))
            //.epsilon(errTol));
        }
        // ("distPerp = 0.1 > D+ell0, single peaked")
        distPerp = 0.06;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
        // ("distPerp = 0.06 < D+ell0, double peaked")
        distPerp = 0.1;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
    }

    SECTION("Test MEDIUM spring") {
        constexpr double errTol = 1e-4;
        const double M = 1.0 / (0.00411);

        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        LookupTable LUT(&lut_filler);
        const double lUB = LUT.getLUCutoff();

        distPerp = 0.04;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
        // ("distPerp = 0.1 > D+ell0, single peaked")
        distPerp = 0.06;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
        // ("distPerp = 0.06 < D+ell0, double peaked")
        distPerp = 0.1;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
    }

    SECTION("Test STIFF spring") {
        constexpr double errTol = 1e-4;
        const double M = 10 / (0.00411);

        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
        LookupTable LUT(&lut_filler);
        const double lUB = LUT.getLUCutoff();

        distPerp = 0.04;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
        // ("distPerp = 0.1 > D+ell0, single peaked")
        distPerp = 0.06;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
        // ("distPerp = 0.06 < D+ell0, double peaked")
        distPerp = 0.1;
        for (double fact = 0.1; fact <= 1; fact += 0.1) {
            CHECK(LUT.Lookup(distPerp, fact * lUB) ==
                  Approx(fdep_integral(distPerp, 0, fact * lUB, M, e_fact,
                                       fdep_length, freelength))
                      .epsilon(errTol));
        }
    }
}

/*
 *TEST_CASE("Lookup table test manual medium spring REL error", "[lookup]")
 *{
 *    // integrated by mathematica
 *    LookupTable LUT(&lut_filler);
 *    const double D = 0.024;
 *    constexpr double errTol = RELTOL;
 *
 *    double distPerp = 0;
 *    lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
 *
 *    distPerp = 0.2;
 *    // ("distPerp = 0.2 > D+ell0, single peaked")
 *    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.0722077) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.142839) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 0.210412) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 0.273623) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 0.383039) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 0.466375) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 0.523889) <
 *errTol); CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 0.55967) <
 *errTol);
 *
 *    distPerp = 0.08;
 *    // "distPerp = 0.08 > D+ell0, single peaked"/D,
 *    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.497588) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.993608) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 1.48554) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 1.96929) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 2.88784) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 3.6925) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 4.3332) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 4.78986) < errTol);
 *
 *    distPerp = 0.06;
 *    // "distPerp = 0.06 < D+ell0, double peaked"/D,
 *    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.488864) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.981139) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 1.47815) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 1.97788) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 2.96052) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 3.85857) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 4.59864) < errTol);
 *    CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 5.14058) < errTol);
 *}
 */

/*
 *TEST_CASE("REVERSE Lookup table test manual medium spring REL error",
 *          "[REVERSE lookup]") {
 *    // integrated by mathematica
 *    LookupTable LUT(&lut_filler);
 *    const double D = 0.024;
 *
 *    double distPerp = 0;
 *    lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
 *
 *    double tol = RELTOL * REVERSEFAC;
 *
 *    distPerp = 0.1;
 *    // ("distPerp = 0.1 > D+ell0, single peaked")
 *    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0) / D, 0.0) < tol);
 *    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.0460519) / D, 0.05) <
 *tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.322128) / D, 0.35)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.459824) / D, 0.5)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.915356) / D, 1.0)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 1.36196) / D, 1.5)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 1.79446) / D, 2.0)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 2.20718) / D, 2.5)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 3.27015) / D, 4.0)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 3.79115) / D, 5.0)
 *< tol); CHECK(relError(LUT.ReverseLookup(distPerp, D * 4.15242) / D, 6.0)
 *< tol);
 *    // CHECK(relError(LUT.ReverseLookup(distPerp / D, 4.37561), 7.0) <
 *tol);
 *}
 */

/*
 *TEST_CASE("REVERSE Lookup table test", "[REVERSE lookup]") {
 *    const double D = 0.024;
 *    const double e_fact = .5;
 *    const double fdep_length = .01;
 *    const double freelength = 0.05;
 *    const double ell0 = freelength / D;
 *
 *    SECTION("Test SOFT spring") {
 *        const double M = 0.1 / (2 * 0.00411);
 *
 *        LookupTable LUT(&lut_filler);
 *        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
 *
 *        double distPerp = 0;
 *        distPerp = 0.2;
 *        // ("distPerp = 0.2 > D+ell0, single peaked")
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
 *                            REVERSEFAC));
 *        }
 *        // ("distPerp = 0.1 > D+ell0, single peaked")
 *        distPerp = 0.1;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
 *                            REVERSEFAC));
 *        }
 *        // ("distPerp = 0.06 < D+ell0, double peaked")
 *        distPerp = 0.06;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
 *                            REVERSEFAC));
 *        }
 *    }
 *
 *    SECTION("Test MEDIUM spring") {
 *        const double M = 1.0 / (2 * 0.00411);
 *
 *        LookupTable LUT(&lut_filler);
 *        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
 *
 *        double distPerp = 0;
 *        distPerp = 0.2;
 *        // ("distPerp = 0.2 > D+ell0, single peaked")
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 3;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound *
 *D));
 *        }
 *        // ("distPerp = 0.1 > D+ell0, single peaked")
 *        distPerp = 0.1;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound *
 *D));
 *        }
 *        // ("distPerp = 0.06 < D+ell0, double peaked")
 *        distPerp = 0.06;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.2) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound *
 *D));
 *        }
 *    }
 *
 *    SECTION("Test STIFF spring") {
 *        const double M = 10 / (2 * 0.00411);
 *
 *        LookupTable LUT(&lut_filler);
 *        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
 *
 *        double distPerp = 0;
 *        distPerp = 0.1;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.1) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
 *                            REVERSEFAC));
 *        }
 *        // ("distPerp = 0.06 < D+ell0, double peaked")
 *        distPerp = 0.06;
 *        for (double sbound = 0; sbound < LUT.getNonDsbound() / 2;
 *             sbound += 0.1) {
 *            double val = fdep_integral(distPerp / D, 0, sbound / D, e_fact,
 *                                       fdep_length / D, M * D * D, ell0);
 *            CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
 *                            REVERSEFAC));
 *        }
 *    }
 *}
 */

/*
 *TEST_CASE("REVERSE binary Lookup test", "[REVERSE binary]") {
 *
 *    const double D = 0.024;
 *    const double e_fact = .5;
 *    const double fdep_length = .01;
 *    const double freelength = 0.05;
 *    const double ell0 = freelength / D;
 *
 *    LookupTable LUT(&lut_filler);
 *
 *    double rowIndexMax = LUT.distPerpGridNumber;
 *    double colIndexMax = LUT.sboundGridNumber;
 *    double distPerpSpacing;
 *
 *    SECTION("Soft spring") {
 *        const double M = .1 / (2 * 0.00411);
 *        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
 *
 *        distPerpSpacing = LUT.distPerpGridSpacing;
 *        for (int i = 0; i < rowIndexMax - 2; ++i) {
 *            double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
 *            double C1 =
 *                LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
 *            double Cavg = .5 * (C0 + C1);
 *            double distPerpAvg = distPerpSpacing * (i + .5) * D;
 *            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
 *            double Cfdep_integral =
 *                fdep_integral(distPerpAvg / D, 0, sbound, M * D * D, ell0) *
 *D; CHECK(errorPass(Cfdep_integral, Cavg));
 *        }
 *    }
 *
 *    SECTION("Medium spring") {
 *        const double M = .1 / (2 * 0.00411);
 *        lut_filler.Init(M, e_fact, fdep_length, freelength, D);
 *
 *        distPerpSpacing = LUT.distPerpGridSpacing;
 *        distPerpSpacing = LUT.distPerpGridSpacing;
 *        for (int i = 0; i < rowIndexMax - 2; ++i) {
 *            double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
 *            double C1 =
 *                LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
 *            double Cavg = .5 * (C0 + C1);
 *            double distPerpAvg = distPerpSpacing * (i + .5) * D;
 *            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
 *            double Cfdep_integral =
 *                fdep_integral(distPerpAvg / D, 0, sbound, M * D * D, ell0) *
 *D; CHECK(errorPass(Cfdep_integral, Cavg));
 *        }
 *    }
 *
 *    SECTION("Stiff spring") {
 *        alpha = 10. / (2 * 0.00411);
 *        M = alpha * D * D;
 *        lut_filler.Init(alpha, freelength, D);
 *        distPerpSpacing = LUT.distPerpGridSpacing;
 *        for (int i = 0; i < rowIndexMax - 2; ++i) {
 *            double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
 *            double C1 =
 *                LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
 *            double Cavg = .5 * (C0 + C1);
 *            double distPerpAvg = distPerpSpacing * (i + .5) * D;
 *            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
 *            double Cfdep_integral =
 *                fdep_integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
 *            CHECK(errorPass(Cfdep_integral, Cavg));
 *        }
 *    }
 *}
 *
 *TEST_CASE("Test the calculation of binding volume.", "[bind volume]") {
 *    const double D = 0.024;
 *    const double freelength = 0.05;
 *    const double ell0 = freelength;
 *    double alpha = 1. / (2 * 0.00411);
 *
 *    for (double i = 0.1; i < 1.0; i += .1) {
 *        LookupTable LUT(&lut_filler);
 *        lut_filler.Init(alpha * i, freelength, D);
 *        LUT.calcBindVol();
 *        double bind_vol =
 *            fdep_bind_vol_integral(0, LUT.getLUCutoff(), i * alpha, ell0);
 *        REQUIRE(LUT.getBindVolume() == Approx(bind_vol).epsilon(1e-8));
 *    }
 *}
 */
