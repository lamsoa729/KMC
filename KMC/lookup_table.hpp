/**
 * @file lookup_table.hpp
 * @author Wen Yan
 * @brief 2D lookup table for probability CDF calculations
 * @version 0.1
 * @date 2019-03-08
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef LOOKUP_TABLE_HPP_
#define LOOKUP_TABLE_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "integrals.hpp"
#include "macros.hpp"
/**
 * @brief 2D Lookup Table for the dimensionless functions
 *
 * Init() is not thread-safe
 * Lookup() and ReverseLookup() are thread safe
 *
 * Both interface functions and internal caculations are dimensionless
 */
class LookupTable {
  private:
    double lUB_; ///< Upper bound distance of lookup table, dimensionless
    double M_;   ///< (1-\lambda)\kappa\beta/2 * D^2, dimensionless
    double ell0; ///< \ell_0/D, dimensionless
    double bind_vol_ = 1; ///< Volume that head can bind within, dimensionless

  public:
    static constexpr double small = 1e-4;
    double D_; ///< dimensional rod Diameter, length scale

    int distPerpGridNumber;
    double distPerpGridSpacing;       ///< dimensionless
    double distPerpGridSpacingInv;    ///< dimensionless
    std::vector<double> distPerpGrid; ///< dimensionless vertical direction

    int sboundGridNumber;
    double sboundGridSpacing;       ///< dimensionless
    double sboundGridSpacingInv;    ///< dimensionless
    std::vector<double> sboundGrid; ///< dimensionless horizontal direction

    std::vector<double> table; ///< the 2D matrix of dimensionless values

  public:
    LookupTable() = default;
    ~LookupTable() = default;

    /**
     * @brief Get the rodDiameter length scale
     *
     * @return double
     */
    double getRodDiameter() const { return D_; }

    /**
     * @brief Get the dimensional biniding volume of attaching end
     *
     * @return double
     */
    double getBindVolume() const { return CUBE(D_) * bind_vol_; }

    /**
     * @brief Get the dimensional cutoff distance of the lookup table.
     *
     * @return double
     */
    double getLUCutoff() const { return D_ * lUB_; }

    /**
     * @brief Get the dimensionless sbound for tabulation
     *
     * @return double
     */
    double getNonDsbound() const { return sboundGrid.back(); }
    /**
     * @brief Get the dimensionful sbound for tabulation
     *
     * @return double
     */
    double getDsbound() const { return D_ * sboundGrid.back(); }

    /*! \brief Initialize lookup table
     *
     * Stores precalculated values of the integral
     * \int_0^y exp{-M_ (sqrt[x^2 + l_min^2] - l_o - D)^2}dx
     * To do this, grid spacing are calculated based on the value of M_ and
     * freeLength. The upper bound lUB is calculated so that the maximum change
     * at the end of the integral is less than 10^{-4}. Everything is
     * non-dimensionalized by rodD.
     *
     * \param M Exponential prefactor with the dimensions of L^{-2}.
     * \param freeLength Rest length of binding crosslink.
     * \param rodD Diameter of rod crosslink is binding to.
     *
     * \return void
     */
    void Init(double M, double freeLength, double rodD) {
        // constexpr double small = 1e-4;
        // setup length scale and dimensionless lengths
        D_ = rodD;
        ell0 = freeLength / rodD;
        M_ = M * rodD * rodD;

        // truncate the integration when integrand < SMALL
        // interation table in dimensionless lengths
        // lUB = sqrt(lm^2 + s^2), dimensionless scaled by D
        lUB_ = sqrt(-log(small) / M_) + ell0;
        // Get binding volume for simulation
        // bind_vol_ = bind_vol_integral(0, lUB_, M_, ell0);
        // printf("lUB_ = %f\n", lUB_);
        // printf("bind_vol_ = %f\n", bind_vol_);

        // step 1 determine grid
        const double distPerpLB = 0;
        const double distPerpUB = lUB_;
        distPerpGridNumber = 256; // grid in dperp
        distPerpGridSpacing =
            (distPerpUB - distPerpLB) / (distPerpGridNumber - 1);

        const double sboundLB = 0;
        const double sboundUB = sqrt(lUB_ * lUB_ - distPerpLB * distPerpLB);
        sboundGridNumber = 256; // grid in s bound
        sboundGridSpacing = (sboundUB - sboundLB) / (sboundGridNumber - 1);

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

    /*! \brief Calculate the binding volume of second head of protein.
     * Will remain 1 otherwise.
     *
     * \return void
     */
    void calcBindVol() { bind_vol_ = bind_vol_integral(0, lUB_, M_, ell0); }

    /*! \brief Set the binding volume of second head of protein.
     * Will remain 1 otherwise.
     *
     * \return void
     */
    void setBindVol(double bind_vol) { bind_vol_ = bind_vol; }

    /*********************
     *  output LookupTable
     *********************/

    /*! \brief Print Table to std::cout
     * Currently 2D only
     *
     * \return void
     */
    void PrintTable() const {
        const int distPerpGridNumber = distPerpGrid.size();
        const int sboundGridNumber = sboundGrid.size();

        // IMPORTANT: Table stores dimensionless val
        for (int i = 0; i < distPerpGridNumber; i++) {
            for (int j = 0; j < sboundGridNumber; j++) {
                std::cout << "(distPerp, sbound) = (" << distPerpGrid[i] << ","
                          << sboundGrid[j]
                          << ") = " << table[i * sboundGridNumber + j]
                          << std::endl;
            }
        }
    }

    /*********************
     *  Look Up
     *********************/

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
    double Lookup(double distPerp, double sbound) const {
        // Non-dimensionalize parameters
        distPerp /= D_;
        sbound /= D_;
        double rowFrac = 0, colFrac = 0;
        int rowIndex = getRowIndex(distPerp, rowFrac);
        int colIndex = getColIndex(sbound, colFrac);

        // clamp to grid bound
        if (rowIndex < 0) {
            fprintf(stderr, "Error: distPerp must be positive value: %g \n",
                    distPerp);
            throw std::runtime_error("distPerp must be positive value");
        }
        if (colIndex < 0) {
            fprintf(stderr, "Error: sbound too small: %g \n", sbound);
            throw std::runtime_error("sbound too small");
        }
        if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
            printf("Warning: distPerp %g very large, clamp to grid UB\n",
                   distPerp);
#endif
            rowIndex = distPerpGridNumber - 2;
            rowFrac = 1;
        }
        if (colIndex > sboundGridNumber - 2) {
            //#ifndef NDEBUG
            //            printf("Warning: sbound %g very large, clamp to grid
            //            UB\n", sbound);
            //#endif
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
                     + table[getTableIndex(rowIndex + 1, colIndex + 1)] *
                           colFrac * rowFrac; //

        return val * D_; // IMPORTANT: re-dimensionalize value
    }

    /******************
     * Invert Lookup
     ******************/

    /*! \brief Invert the integral represented by lookup table so given a
     * distance away and value of integral you retrieve the upper limit of
     * integration
     *
     * \param distPerp l_min in the above integral. The distance away from the
     * infinite carrier line.
     * \param val Resultant value of integral
     *
     * \return double Dimensional upper bound of integral
     */
    double ReverseLookup(double distPerp, double val) const {
        if (val == 0) {
            return 0;
        } else { // Non-dimensionalize inputs
            distPerp /= D_;
            val /= D_;
        }

        double sbound;
        double rowFrac = 0;
        int rowIndex = getRowIndex(distPerp, rowFrac);

        // distPerp check
        if (rowIndex < 0) {
            fprintf(stderr, "Error: distPerp must be positive value: %g \n",
                    distPerp);
            throw std::runtime_error("distPerp must be positive value.");
        }
        if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
            printf(
                "Warning: distPerp %g very large, clamp to grid upper limit\n",
                distPerp);
#endif
            rowIndex = distPerpGridNumber - 2;
            rowFrac = 1;
        }

        // Every row of table starts from 0, for sbound=0
        // val should be positive
        if (val < 0) {
#ifndef NDEBUG
            fprintf(stderr, "Error: val should be positive: %g \n", val);
#endif
            throw std::runtime_error("val should be positive.");
        }

        // reverse lookup val on row rowIndex

        // Needed only for quadratic interpolation
        // int index0LB;
        // int index0UB;
        // if (rowIndex > 0) {
        //    index0LB = getTableIndex(rowIndex - 1, 0);
        //    index0UB = getTableIndex(rowIndex - 1, sboundGridNumber - 1);
        //} else {
        //    index0LB = getTableIndex(rowIndex, 0);
        //    index0UB = getTableIndex(rowIndex, sboundGridNumber - 1);
        //}

        const int index1LB = getTableIndex(rowIndex, 0);
        const int index1UB = getTableIndex(rowIndex, sboundGridNumber - 1);
        const int index2LB = getTableIndex(rowIndex + 1, 0);
        const int index2UB = getTableIndex(rowIndex + 1, sboundGridNumber - 1);

        // reverse lookup
        // auto lower0 = std::lower_bound(table.begin() + index0LB,
        // table.begin() + index0UB, val);
        // reverse lookup val on row rowIndex
        auto lower1 = std::lower_bound(table.begin() + index1LB,
                                       table.begin() + index1UB, val);
        // reverse lookup val on row rowIndex+1
        auto lower2 = std::lower_bound(table.begin() + index2LB,
                                       table.begin() + index2UB, val);

        // find cross point of distPerp and two points
        // assert(lower0 - table.begin() >= index0LB);
        assert(lower1 - table.begin() >= index1LB);
        assert(lower2 - table.begin() >= index2LB);

        /**
         * value on table grids
         * rowIndex   --------+--------+---
         *             colIndexm colIndexm+1
         *
         * rowIndex+1 -----------+--------+
         *                  colIndexp colIndexp+1
         *
         * Off set because different distPerps have different cutoff values
         */
        // const int colIndex0 = lower0 - 1 - table.begin() - index1LB;
        const int colIndexm = lower1 - 1 - table.begin() - index1LB;
        const int colIndexp = lower2 - 1 - table.begin() - index2LB;
        // double sbound0;
        double sboundm, sboundp;
        int out_of_range_case = 0;

        //        if (lower0 == table.begin() + index1UB) {
        //#ifndef NDEBUG
        //            printf("Warning: val %g too large for row0 lookup with max
        //            of %g. "
        //                   "Setting to a max sbound of %g \n",
        //                   val, *(lower0), sboundGrid[colIndex0]);
        //#endif
        //            sbound0 = sboundGrid[colIndex0];
        //        } else {
        //            double val0A = *(lower0 - 1);
        //            double val0B = *(lower0);
        //            sbound0 = sboundGrid[colIndex0] +
        //                      (val - val0A) / (val0B - val0A) *
        //                      sboundGridSpacing;
        //        }
        // Interpolate first row in the sbound direction
        if (lower1 == table.begin() + index1UB) {
#ifndef NDEBUG
            printf("Warning: val %g too large for row1 lookup with max of %g. "
                   "Setting sbound max to %g \n",
                   val, *(lower1), sboundGrid[colIndexm]);
#endif
            out_of_range_case = 1;
            sboundm = sboundGrid[colIndexm]; // Go to end of sbound grid
        } else {
            double valmA = *(lower1 - 1);
            double valmB = *(lower1);
            sboundm = sboundGrid[colIndexm] +
                      ((val - valmA) / (valmB - valmA)) * sboundGridSpacing;
        }
        // Interpolate second row in the sbound direction
        if (lower2 == table.begin() + index2UB) {
#ifndef NDEBUG
            printf("Warning: val %g too large for row2 lookup with max of %g. "
                   "Setting sbound max to %g \n",
                   val, *(lower2), sboundGrid[colIndexp]);
#endif
            assert(out_of_range_case != 1); // Something has gone horribly wrong
            out_of_range_case = 2;
            sboundp = sboundGrid[colIndexp]; // Go to end of sbound grid
            // return D * sboundm;
        } else {
            double valpA = *(lower2 - 1);
            double valpB = *(lower2);
            sboundp = sboundGrid[colIndexp] +
                      ((val - valpA) / (valpB - valpA)) * sboundGridSpacing;
        }
        // printf("sboundm = %f\n", sboundm);
        // printf("sboundp = %f\n", sboundp);

        // Interpolate two sbound values in the distPerp direction
        // if (rowIndex == 0) { // Linear interpolation
        switch (out_of_range_case) {
        case 0:
            sbound = sboundm * (1 - rowFrac) + sboundp * rowFrac;
            break;
        case 1:
            sbound =
                ReverseBinaryLookup(rowIndex, rowFrac, sboundp, sboundm, val);
            assert(sbound >= (sboundp) && sbound <= (sboundm));
            break;
        case 2:
            sbound =
                ReverseBinaryLookup(rowIndex, rowFrac, sboundm, sboundp, val);
            assert(sbound <= (sboundp) && sbound >= (sboundm));
            break;
        }
        //} else { // Quadratic interpolation. TODO Needs testing
        //    printf("In quadratic \n");
        //    sbound = .5 * rowFrac *
        //                 ((rowFrac - 1.) * sbound0 + (rowFrac + 1.) *
        //                 sboundp) +
        //             (1 - rowFrac) * (1 + rowFrac) * sboundm;
        //}

#ifndef NDEBUG
        double UB_val1 = *(table.begin() + index1UB);
        double UB_val2 = *(table.begin() + index2UB);
        printf("Reverse lookup input = %g, Upper limits = %g or %g \n", val,
               UB_val1, UB_val2);
        printf("Reverse lookup output = %g, expected range = (%g, %g) \n",
               sbound, sboundm, sboundp);
#endif
        return D_ * sbound;
    }

    double ReverseBinaryLookup(const int rowMin, const double rowFrac,
                               const double sMin, const double sMax,
                               const double val) const {
        // Isn't used, just needed for the getColIndex function
        double colFrac = 0;

        int colMin = getColIndex(sMin, colFrac);
        int colMax = getColIndex(sMax, colFrac);
        // printf("colMin init = %d\n", colMin);
        // printf("colMax init= %d\n", colMax);
        assert(colMin >= 0);
        assert(colMax < sboundGridNumber);
        // Simple binary search using linear interpolation to find value
        while ((colMax - 1) > colMin) {
            double avg = (colMax + colMin) * .5;
            int colAvg = floor(avg);
            double colVal =
                table[getTableIndex(rowMin, colAvg)] * (1 - rowFrac)  //
                + table[getTableIndex(rowMin + 1, colAvg)] * rowFrac; //
            // printf("colAvg = %d, colVal = %f, val = %f\n", colAvg, colVal,
            // val);

            // Accuracy of table is set to 1e-4, use 1e-5 to set colMax value.
            // Otherwise, you will always reach the end of column
            if ((val - colVal) < 1e-5)
                colMax = colAvg;
            else if (colVal < val)
                colMin = colAvg;
        }
        // printf("colMin final = %d\n", colMin);
        // printf("colMax final = %d\n", colMax);
        // Interpolate using function of the line
        double valMin = table[getTableIndex(rowMin, colMin)] * (1 - rowFrac)  //
                        + table[getTableIndex(rowMin + 1, colMin)] * rowFrac; //
        double valMax = table[getTableIndex(rowMin, colMax)] * (1 - rowFrac)  //
                        + table[getTableIndex(rowMin + 1, colMax)] * rowFrac; //

        // Linear interpolation with known values
        double sbound = (val - valMin) / (valMax - valMin) * sboundGridSpacing +
                        sboundGrid[colMin];
        // Make sure sbound is not calculated to be outside range [sMin, sMax]
        if (sbound > sMax)
            return sMax;
        else if (sbound < sMin)
            return sMin;
        else
            return sbound;
    }

    inline int getTableIndex(int row, int col) const {
        assert(row < distPerpGridNumber && row >= 0);
        assert(col < sboundGridNumber && col >= 0);
        return row * sboundGridNumber + col;
    }

  private:
    /*! \brief Fills the lookup table with values after grid has been
     * created.
     *
     * Table is a flattened 2D array with first index corresponding to
     * different distances head is away from rod and second index is defines
     * the bounds of where that crosslinking head can bind.
     *
     * \return void
     */
    void FillMatrix() {
        table.resize(distPerpGridNumber * sboundGridNumber, 0);

        //// boost integration parameters
        // const int max_depth = 10; // maximum number of interval splittings
        // const double tol = 1e-8;  // maximum relative error

        // row major
        // i is slow changing, should be distPerpGrid
        // j is fast changing, should be sboundGrid

        // IMPORTANT: Table stores dimensionless val
        for (int i = 0; i < distPerpGridNumber; i++) {
            for (int j = 0; j < sboundGridNumber; j++) {
                const double sbound = sboundGrid[j];
                assert(!(sbound < 0));
                double result = integral(distPerpGrid[i], 0, sbound, M_, ell0);
                table[i * sboundGridNumber + j] = result;
            }
        }
    }
    /**
     * @brief get row index with distPerp
     *
     * @param distPerp dimensionless
     * @return int
     */
    inline int getRowIndex(double distPerp, double &rowFrac) const {
        const double x = (distPerp - distPerpGrid[0]) * distPerpGridSpacingInv;
        int index = floor(x);
        rowFrac = x - index;
        return index;
    }

    /**
     * @brief Get col index with sbound
     *
     * @param sbound dimensionless
     * @return int
     */
    inline int getColIndex(double sbound, double &colFrac) const {
        const double x = (sbound - sboundGrid[0]) * sboundGridSpacingInv;
        int index = floor(x);
        colFrac = x - index;
        return index;
    }
};

#endif
