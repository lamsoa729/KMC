#include "lookup_table.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>

#include <boost/math/quadrature/gauss_kronrod.hpp>

/*! \brief Initialize lookup table
 *
 * Stores precalculated values of the integral
 * \int_0^y exp{-M_ (sqrt[x^2 + l_min^2] - l_o - D)^2}dx
 * To do this, grid spacing are calculated based on the value of M_ and
 * freeLength. The upper bound lUB is calculated so that the maximum change at
 * the end of the integral is less than 10^{-5}. Everything is
 * non-dimensionalized by rodD.
 *
 * \param M_ Exponential prefactor with the dimensions of L^{-2}.
 * \param freeLength Rest length of binding crosslink.
 * \param rodD Diameter of rod crosslink is binding to.
 *
 * \return void
 */
void LookupTable::Init(double M_, double freeLength, double rodD) {
    // setup length scale and dimensionless lengths
    D = rodD;
    ell0 = freeLength / rodD;
    M = M_ * rodD * rodD;

    // truncate the integration when integrand < SMALL
    // interation table in dimensionless lengths
    // lUB = sqrt(lm^2 + s^2), dimensionless scaled by D
    constexpr double SMALL = 1e-5;
    const double lUB = sqrt(-log(SMALL) / M) + 1 + ell0;

    // step 1 determine grid
    const double distPerpLB =
        0; // allow min distPerp to be half D, some overlap
    const double distPerpUB = lUB;
    distPerpGridNumber = 48; // grid in dperp
    double distPerpGridSpacing =
        (distPerpUB - distPerpLB) / (distPerpGridNumber - 1);

    const double sboundLB = 0;
    const double sboundUB = sqrt(lUB * lUB - distPerpLB * distPerpLB);
    sboundGridNumber = 64; // grid in s bound
    double sboundGridSpacing = (sboundUB - sboundLB) / (sboundGridNumber - 1);

    // step 2 init grid
    distPerpGrid.resize(distPerpGridNumber);
    for (int i = 0; i < distPerpGridNumber; i++) {
        distPerpGrid[i] = i * distPerpGridSpacing + distPerpLB;
    }
    sboundGrid.resize(sboundGridNumber);
    for (int i = 0; i < sboundGridNumber; i++) {
        sboundGrid[i] = i * sboundGridSpacing + sboundLB;
    }
    distPerpGridSpacingInv = 1. / distPerpGridSpacing;
    sboundGridSpacingInv = 1. / sboundGridSpacing;

    // step 3 tabulate the matrix
    FillMatrix();
}

/*! \brief Fills the lookup table with values after grid has been created.
 *
 * Table is a flattened 2D array with first index corresponding to different
 * distances head is away from rod and second index is defines the bounds of
 * where that crosslinking head can bind.
 *
 * \return void
 */
void LookupTable::FillMatrix() {
    table.resize(distPerpGridNumber * sboundGridNumber, 0);

    // boost integration parameters
    const int max_depth = 10; // maximum number of interval splittings
    const double tol = 1e-8;  // maximum relative error

    // row major
    // i is slow changing, should be distPerpGrid
    // j is fast changing, should be sboundGrid

    // IMPORTANT: Table stores dimensionless val
    for (int i = 0; i < distPerpGridNumber; i++) {
        for (int j = 0; j < sboundGridNumber; j++) {
            const double sbound = sboundGrid[j];
            assert(!(sbound < 0));
            double result = integral(distPerpGrid[i], 0, sbound, M, ell0);
            table[i * sboundGridNumber + j] = result;
        }
    }
}

/*! \brief Print out table for testing purposes
 *
 * \return void
 */
void LookupTable::PrintTable() const {
    const int distPerpGridNumber = distPerpGrid.size();
    const int sboundGridNumber = sboundGrid.size();

    // IMPORTANT: Table stores dimensionless val
    for (int i = 0; i < distPerpGridNumber; i++) {
        for (int j = 0; j < sboundGridNumber; j++) {
            std::cout << "(distPerp, sbound) = (" << distPerpGrid[i] << ","
                      << sboundGrid[j]
                      << ") = " << table[i * sboundGridNumber + j] << std::endl;
        }
    }
}

/*! \brief Lookup the value of the integral given two parameters.
 *
 * $\int_0^y exp{-M_ (sqrt[x^2 + l_min^2] - l_o - D)^2}dx$
 *
 * \param distPerp l_min in the above integral. The distance away from the
 * infinite carrier line.
 * \param sbound y in the above integral. The limit of
 * integration.
 *
 * \return double Dimensional value of integral calculated
 */
double LookupTable::Lookup(double distPerp, double sbound) const {
    // Non-dimensionalize parameters
    distPerp /= D;
    sbound /= D;
    double rowFrac = 0, colFrac = 0;
    int rowIndex = getRowIndex(distPerp, rowFrac);
    int colIndex = getColIndex(sbound, colFrac);

    // clamp to grid bound
    if (rowIndex < 0) {
        printf("Error: distPerp too small: %g \n", distPerp);
        exit(1);
    }
    if (colIndex < 0) {
        printf("Error: sbound too small: %g \n", sbound);
        exit(1);
    }
    if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: distPerp %g very large, clamp to grid UB\n", distPerp);
#endif
        rowIndex = distPerpGridNumber - 2;
        rowFrac = 1;
    }
    if (colIndex > sboundGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: sbound %g very large, clamp to grid UB\n", sbound);
#endif
        colIndex = sboundGridNumber - 2;
        colFrac = 1;
    }

    // linear interpolation using surrounding grid values
    double val = table[getTableIndex(rowIndex, colIndex)] * (1 - colFrac) *
                     (1 - rowFrac) //
                 + table[getTableIndex(rowIndex, colIndex + 1)] * colFrac *
                       (1 - rowFrac) //
                 + table[getTableIndex(rowIndex + 1, colIndex)] *
                       (1 - colFrac) * rowFrac //
                 + table[getTableIndex(rowIndex + 1, colIndex + 1)] * colFrac *
                       rowFrac; //

    return val * D; // IMPORTANT: re-dimensionalize value
}

/*! \brief Invert the integral represented by lookup table so given a distance
 * away and value of integral you retrieve the upper limit of integration
 *
 * \param distPerp l_min in the above integral. The distance away from the
 * infinite carrier line.
 * \param val Resultant value of integral
 *
 * \return double Dimensional upper bound of integral
 */
double LookupTable::ReverseLookup(double distPerp, double val) const {
    if (val == 0) {
        return 0;
    } else { // Non-dimensionalize inputs
        distPerp /= D;
        val /= D;
    }

    double sbound;
    double rowFrac = 0;
    int rowIndex = getRowIndex(distPerp, rowFrac);

    // distPerp check
    if (rowIndex < 0) {
        printf("Error: distPerp too small: %g \n", distPerp);
        exit(1);
    }
    if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: distPerp %g very large, clamp to grid UB\n", distPerp);
#endif
        rowIndex = distPerpGridNumber - 2;
        rowFrac = 1;
    }

    // Every row of table starts from 0, for sbound=0
    // val should be positive
    if (val < 0) {
#ifndef NDEBUG
        printf("Error: val should be positive: %g \n", val);
#endif
        exit(1);
    }
    double colFrac = 0;

    // reverse lookup val on row rowIndex
    const int index1LB = getTableIndex(rowIndex, 0);
    const int index1UB = getTableIndex(rowIndex, sboundGridNumber - 1);
    const int index2LB = getTableIndex(rowIndex + 1, 0);
    const int index2UB = getTableIndex(rowIndex + 1, sboundGridNumber - 1);

    auto lower1 = std::lower_bound(table.begin() + index1LB,
                                   table.begin() + index1UB, val);
    // reverse lookup val on row rowIndex+1
    auto lower2 = std::lower_bound(table.begin() + index2LB,
                                   table.begin() + index2UB, val);

    if (lower1 == table.begin() + index1UB ||
        lower2 == table.begin() + index2UB) {
        sbound = sboundGrid.back();
#ifndef NDEBUG
        printf("Warning: val %g too large, setting sbound to max %g\n", val,
               sbound);
#endif
        return D * sbound; // Re-dimensionalize
    }

    // find cross point of distPerp and two points
    assert(lower1 - table.begin() >= index1LB);
    assert(lower2 - table.begin() >= index2LB);

    /**
     * value on table grids
     * rowIndex --------+--------+---
     *             colIndexm colIndexm+1
     *
     * rowIndex+1 -----------+--------+
     *             colIndexp colIndexp+1
     */
    const int colIndexm = lower1 - 1 - table.begin() - index1LB;
    const int colIndexp = lower2 - 1 - table.begin() - index2LB;
    double valmA = *(lower1 - 1);
    double valmB = *(lower1);
    double valpA = *(lower2 - 1);
    double valpB = *(lower2);
    double sboundGridSpacing = sboundGrid[1] - sboundGrid[0];
    double sboundm = sboundGrid[colIndexm] +
                     (val - valmA) / (valmB - valmA) * sboundGridSpacing;
    double sboundp = sboundGrid[colIndexp] +
                     (val - valpA) / (valpB - valpA) * sboundGridSpacing;
    sbound = sboundm * (1 - rowFrac) + sboundp * rowFrac;
    return D * sbound;
}
