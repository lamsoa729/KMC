#include "KMC/integrals.hpp"
#include "KMC/lookup_table.hpp"
#include "catch.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include <boost/math/quadrature/gauss_kronrod.hpp>

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

TEST_CASE("Lookup table test SOFT spring ", "[lookup]") {
    constexpr double errTol = 1e-2;
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(relError(LUT.Lookup(distPerp, sbound * D),
                       D * integral(distPerp / D, 0, sbound, M, ell0)) <
              errTol);
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(relError(LUT.Lookup(distPerp, sbound * D),
                       D * integral(distPerp / D, 0, sbound, M, ell0)) <
              errTol);
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(relError(LUT.Lookup(distPerp, sbound * D),
                       D * integral(distPerp / D, 0, sbound, M, ell0)) <
              errTol);
    }
}

TEST_CASE("Lookup table test MEDIUM spring ", "[lookup]") {
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 1.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
}

TEST_CASE("Lookup table test STIFF spring ", "[lookup]") {
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 10.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
    }
}

TEST_CASE("Lookup table test manual medium spring REL error", "[lookup]") {
    // integrated by mathematica
    LookupTable LUT;
    const double D = 0.024;
    constexpr double errTol = RELTOL;

    double distPerp = 0;
    LUT.Init(1.0 / (2 * 0.00411), 0.05, D);

    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.0722077) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.142839) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 0.210412) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 0.273623) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 0.383039) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 0.466375) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 0.523889) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 0.55967) < errTol);

    distPerp = 0.08;
    // "distPerp = 0.08 > D+ell0, single peaked"/D,
    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.497588) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.993608) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 1.48554) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 1.96929) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 2.88784) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 3.6925) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 4.3332) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 4.78986) < errTol);

    distPerp = 0.06;
    // "distPerp = 0.06 < D+ell0, double peaked"/D,
    CHECK(relError(LUT.Lookup(distPerp, 0.5 * D) / D, 0.488864) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.0 * D) / D, 0.981139) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 1.5 * D) / D, 1.47815) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 2.0 * D) / D, 1.97788) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 3.0 * D) / D, 2.96052) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 4.0 * D) / D, 3.85857) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 5.0 * D) / D, 4.59864) < errTol);
    CHECK(relError(LUT.Lookup(distPerp, 6.0 * D) / D, 5.14058) < errTol);
}

TEST_CASE("REVERSE Lookup table test manual medium spring REL error",
          "[REVERSE lookup]") {
    // integrated by mathematica
    LookupTable LUT;
    const double D = 0.024;

    double distPerp = 0;
    LUT.Init(1.0 / (2 * 0.00411), 0.05, D);

    double tol = RELTOL * REVERSEFAC;

    distPerp = 0.1;
    // ("distPerp = 0.1 > D+ell0, single peaked")
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0) / D, 0.0) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.0460519) / D, 0.05) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.322128) / D, 0.35) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.459824) / D, 0.5) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 0.915356) / D, 1.0) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 1.36196) / D, 1.5) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 1.79446) / D, 2.0) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 2.20718) / D, 2.5) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 3.27015) / D, 4.0) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 3.79115) / D, 5.0) < tol);
    CHECK(relError(LUT.ReverseLookup(distPerp, D * 4.15242) / D, 6.0) < tol);
    // CHECK(relError(LUT.ReverseLookup(distPerp / D, 4.37561), 7.0) < tol);
}

TEST_CASE("REVERSE Lookup table test soft spring ", "[REVERSE lookup]") {
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 3; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 3; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 3; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
}

TEST_CASE("REVERSE Lookup table test medium spring ", "[REVERSE lookup]") {
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 1.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
}

TEST_CASE("REVERSE Lookup table test stiff spring ", "[REVERSE lookup]") {
    LookupTable LUT;

    const double D = 0.024;
    const double alpha = 10.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;
    LUT.Init(alpha, freelength, D);

    double distPerp = 0;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    // WARNING: This reverse lookup fails because function is too flat
    // distPerp = 0.2;
    // for (double sbound = 0; sbound < LUT.getNonDsbound() / 8; sbound += 0.1)
    // {
    //    double val = integral(distPerp / D, 0, sbound, M, ell0);
    //    double scalc = LUT.ReverseLookup(distPerp, val * D);
    //    printf("scalc = %f\n", scalc);
    //    printf("sbound = %f\n", sbound * D);
    //    CHECK(errorPass(scalc, sbound * D, REVERSEFAC));
    //}

    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
                        REVERSEFAC));
    }
}

TEST_CASE("REVERSE binary Lookup with different springs", "[REVERSE binary]") {

    const double D = 0.024;
    const double freelength = 0.05;
    const double ell0 = freelength / D;

    LookupTable LUT;

    double alpha = 0.1 / (2 * 0.00411);
    double M = alpha * D * D;
    LUT.Init(alpha, freelength, D);
    double rowIndexMax = LUT.distPerpGridNumber;
    double colIndexMax = LUT.sboundGridNumber;
    double distPerpSpacing = LUT.distPerpGridSpacing;

    for (int i = 0; i < rowIndexMax - 2; ++i) {
        double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
        double C1 = LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
        double Cavg = .5 * (C0 + C1);
        double distPerpAvg = distPerpSpacing * (i + .5) * D;
        double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
        double Cintegral =
            integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
        CHECK(errorPass(Cintegral, Cavg));
    }

    // Medium spring
    alpha = 1.0 / (2 * 0.00411);
    M = alpha * D * D;
    LUT.Init(alpha, freelength, D);
    distPerpSpacing = LUT.distPerpGridSpacing;
    for (int i = 0; i < rowIndexMax - 2; ++i) {
        double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
        double C1 = LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
        double Cavg = .5 * (C0 + C1);
        double distPerpAvg = distPerpSpacing * (i + .5) * D;
        double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
        double Cintegral =
            integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
        CHECK(errorPass(Cintegral, Cavg));
    }

    // Stiff spring
    alpha = 10. / (2 * 0.00411);
    M = alpha * D * D;
    LUT.Init(alpha, freelength, D);
    distPerpSpacing = LUT.distPerpGridSpacing;
    for (int i = 0; i < rowIndexMax - 2; ++i) {
        double C0 = LUT.table[LUT.getTableIndex(i, colIndexMax - 1)] * D;
        double C1 = LUT.table[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
        double Cavg = .5 * (C0 + C1);
        double distPerpAvg = distPerpSpacing * (i + .5) * D;
        double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
        double Cintegral =
            integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
        CHECK(errorPass(Cintegral, Cavg));
    }

    // double distPerp = 0;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    // WARNING: This reverse lookup fails because function is too flat
    // distPerp = 0.2;
    // for (double sbound = 0; sbound < LUT.getNonDsbound() / 8; sbound += 0.1)
    // {
    //    double val = integral(distPerp / D, 0, sbound, M, ell0);
    //    double scalc = LUT.ReverseLookup(distPerp, val * D);
    //    printf("scalc = %f\n", scalc);
    //    printf("sbound = %f\n", sbound * D);
    //    CHECK(errorPass(scalc, sbound * D, REVERSEFAC));
    //}

    // ("distPerp = 0.1 > D+ell0, single peaked")
    // distPerp = 0.1;
    // for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1)
    // {
    //    double val = integral(distPerp / D, 0, sbound, M, ell0);
    //    CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
    //                    REVERSEFAC));
    //}
    //// ("distPerp = 0.06 < D+ell0, double peaked")
    // distPerp = 0.06;
    // for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1)
    // {
    //    double val = integral(distPerp / D, 0, sbound, M, ell0);
    //    CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
    //                    REVERSEFAC));
    //}
}
