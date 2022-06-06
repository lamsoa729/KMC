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
#include "lut_filler.hpp"
#include "lut_filler_asym.hpp"
#include "lut_filler_edep.hpp"
#include "lut_filler_fdep.hpp"
#include "lut_filler_2nd_order.hpp"
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
    // All of the variables from lut_filler
    double exp_fact1_;
    double exp_fact2_;
    double e_fact_;
    double fdep_length_;
    double rest_length_;

  protected:
    double UB_; ///< Upper bound distance of lookup table, dimensionless
    double D_;  ///< dimensional rod Diameter, length scale

    double bind_vol_ = 1.; ///< Volume that head can bind within, dimensionful
  public:
    static constexpr double small_ = 1e-4;

    int perp_grid_num_ = 256;
    double perp_spacing_;           ///< dimensionless
    double perp_spacing_inv_;       ///< dimensionless
    std::vector<double> perp_grid_; ///< dimensionless vertical direction

    int para_grid_num_ = 256;
    double para_spacing_;           ///< dimensionless
    double para_spacing_inv_;       ///< dimensionless
    std::vector<double> para_grid_; ///< dimensionless horizontal direction

    std::vector<double> table_; ///< the 2D matrix of dimensionless values

  public:
    LookupTable() = default;
    ~LookupTable() = default;

    double CalcBoltzmann(double distCent);

    /*! \brief Initialize lookup table
     *
     * Stores precalculated values of the integral from lut_filler
     *
     * \param lut_filler Initialized lookup table filler pointer
     *
     * \return void
     */
    LookupTable(LUTFiller *lut_filler, bool use_bind_vol = false) {
        // setup length scale and dimensionless lengths

        // truncate the integration when integrand < SMALL
        // interation table in dimensionless lengths
        exp_fact1_ = lut_filler->getExpFact1();
        exp_fact2_ = lut_filler->getExpFact2();
        e_fact_ = lut_filler->getEFact();
        fdep_length_ = lut_filler->getFDepLength();
        rest_length_ = lut_filler->getRestLength();

        UB_ = lut_filler->getUpperBound();
        D_ = lut_filler->getLengthScale();

        para_grid_num_ = lut_filler->getDistParaGridNum();
        perp_grid_num_ = lut_filler->getDistPerpGridNum();

        para_spacing_ = lut_filler->getDistParaGridSpacing();
        perp_spacing_ = lut_filler->getDistPerpGridSpacing();
        para_spacing_inv_ = 1. / para_spacing_;
        perp_spacing_inv_ = 1. / perp_spacing_;

        lut_filler->FillDistParaGrid(para_grid_);
        lut_filler->FillDistPerpGrid(perp_grid_);

        lut_filler->FillMatrix(table_);
        // If not using binding volume, binding volume will be set to 1 so
        // it does not effect KMC calculations
        if (use_bind_vol) {
            setBindVol(lut_filler->getBindingVolume());
        } else {
            setDBindVol(1.);
        }
    }

    /**
     * @brief Calculate Boltzmann factor for point-like object
     *
     * @return double Boltzmann factor
     */
    inline double calcBoltzmann(double dist_cent) const {
        /* TODO: Make dimensionful unit test for this <30-03-21, ARL> */
        double r = dist_cent / D_;
        // exp_fact1_ used for compressed spring, exp_fact2_ for stretched
        double exp_fact = r < rest_length_ ? exp_fact1_ : exp_fact2_;
        return exp(-exp_fact * ((1. - e_fact_) * SQR(r - rest_length_) -
                                fdep_length_ * (r - rest_length_)));
    }

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
    double getBindVolume() const { return bind_vol_; }

    /**
     * @brief Get the dimensional cutoff distance of the lookup table.
     *
     * @return double
     */
    double getLUCutoff() const { return D_ * UB_; }

    /**
     * @brief Get the dimensionless sbound for tabulation
     *
     * @return double
     */
    double getNonDsbound() const { return para_grid_.back(); }

    /**
     * @brief Get the dimensionful sbound for tabulation
     *
     * @return double
     */
    double getDsbound() const { return D_ * para_grid_.back(); }

    /*! \brief Set the binding volume of second head of protein making it
     * dimensionful. Will remain 1 otherwise.
     *
     * \return void
     */
    void setBindVol(double bind_vol) {
        assert(bind_vol > 0);
        bind_vol_ = CUBE(D_) * bind_vol;
    }

    /*! \brief Set dimensionful binding volume of second head of protein.
     * Will remain 1 otherwise.
     *
     * \return void
     */

    void setDBindVol(double bind_vol) {
        assert(bind_vol > 0);
        bind_vol_ = bind_vol;
    }

    /*********************
     *  output LookupTable
     *********************/

    /*! \brief Print Table to std::cout
     * Currently 2D only
     *
     * \return void
     */
    void PrintTable() const {
        const int perp_grid_num_ = perp_grid_.size();
        const int para_grid_num_ = para_grid_.size();

        // IMPORTANT: Table stores dimensionless val
        for (int i = 0; i < perp_grid_num_; i++) {
            for (int j = 0; j < para_grid_num_; j++) {
                std::cout << "(dist_perp, dist_para) = (" << perp_grid_[i]
                          << "," << para_grid_[j]
                          << ") = " << table_[i * para_grid_num_ + j]
                          << std::endl;
            }
        }
    }

    /*********************
     *  Look Up
     *********************/

    /*! \brief Lookup the value of the integral given two parameters.
     *
     *
     * \param dist_perp The distance away from the
     * infinite carrier line.
     * \param dist_para The limit of integration.
     *
     * \return double Dimensional value of integral calculated
     */
    double Lookup(double dist_perp, double dist_para) const {
        // Non-dimensionalize parameters
        dist_perp /= D_;
        dist_para /= D_;
        double rowFrac = 0, colFrac = 0;
        int rowIndex = getRowIndex(dist_perp, rowFrac);
        int colIndex = getColIndex(dist_para, colFrac);

        // clamp to grid bound
        if (rowIndex < 0) {
            fprintf(stderr, "Error: dist_perp must be positive value: %g \n",
                    dist_perp);
            throw std::runtime_error("dist_perp must be positive value");
        }
        if (colIndex < 0) {
            fprintf(stderr, "Error: sbound too small: %g \n", dist_para);
            throw std::runtime_error("sbound too small");
        }
        if (rowIndex > perp_grid_num_ - 2) {
#ifndef NDEBUG
            printf("Warning: dist_perp %g very large, clamp to grid UB\n",
                   dist_perp);
#endif
            rowIndex = perp_grid_num_ - 2;
            rowFrac = 1;
        }
        if (colIndex > para_grid_num_ - 2) {
            //#ifndef NDEBUG
            //            printf("Warning: sbound %g very large, clamp to
            //            grid UB\n", sbound);
            //#endif
            colIndex = para_grid_num_ - 2;
            colFrac = 1;
        }

        // linear interpolation using surrounding grid values
        double val = table_[getTableIndex(rowIndex, colIndex)] * (1 - colFrac) *
                         (1 - rowFrac) //
                     + table_[getTableIndex(rowIndex, colIndex + 1)] * colFrac *
                           (1 - rowFrac) //
                     + table_[getTableIndex(rowIndex + 1, colIndex)] *
                           (1 - colFrac) * rowFrac //
                     + table_[getTableIndex(rowIndex + 1, colIndex + 1)] *
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
     * \param dist_perp l_min in the above integral. The distance away from
     * the infinite carrier line. \param val Resultant value of integral
     *
     * \return double Dimensional upper bound of integral
     */
    double ReverseLookup(double dist_perp, double val) const {
        if (val == 0) {
            return 0;
        } else { // Non-dimensionalize inputs
            dist_perp /= D_;
            val /= D_;
        }

        double sbound;
        double rowFrac = 0;
        int rowIndex = getRowIndex(dist_perp, rowFrac);

        // dist_perp check
        if (rowIndex < 0) {
            fprintf(stderr, "Error: dist_perp must be positive value: %g \n",
                    dist_perp);
            throw std::runtime_error("dist_perp must be positive value.");
        }
        if (rowIndex > perp_grid_num_ - 2) {
#ifndef NDEBUG
            printf("Warning: dist_perp %g very large, clamp to grid upper "
                   "limit\n",
                   dist_perp);
#endif
            rowIndex = perp_grid_num_ - 2;
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
        //    index0UB = getTableIndex(rowIndex - 1, para_grid_num_ - 1);
        //} else {
        //    index0LB = getTableIndex(rowIndex, 0);
        //    index0UB = getTableIndex(rowIndex, para_grid_num_ - 1);
        //}

        const int index1LB = getTableIndex(rowIndex, 0);
        const int index1UB = getTableIndex(rowIndex, para_grid_num_ - 1);
        const int index2LB = getTableIndex(rowIndex + 1, 0);
        const int index2UB = getTableIndex(rowIndex + 1, para_grid_num_ - 1);

        // reverse lookup
        // auto lower0 = std::lower_bound(table.begin() + index0LB,
        // table.begin() + index0UB, val);
        // reverse lookup val on row rowIndex
        auto lower1 = std::lower_bound(table_.begin() + index1LB,
                                       table_.begin() + index1UB, val);
        // reverse lookup val on row rowIndex+1
        auto lower2 = std::lower_bound(table_.begin() + index2LB,
                                       table_.begin() + index2UB, val);

        // find cross point of dist_perp and two points
        // assert(lower0 - table.begin() >= index0LB);
        assert(lower1 - table_.begin() >= index1LB);
        assert(lower2 - table_.begin() >= index2LB);

        /**
         * value on table grids
         * rowIndex   --------+--------+---
         *             colIndexm colIndexm+1
         *
         * rowIndex+1 -----------+--------+
         *                  colIndexp colIndexp+1
         *
         * Off set because different dist_perps have different cutoff values
         */
        // const int colIndex0 = lower0 - 1 - table.begin() - index1LB;
        const int colIndexm = lower1 - 1 - table_.begin() - index1LB;
        const int colIndexp = lower2 - 1 - table_.begin() - index2LB;
        // double sbound0;
        double sboundm, sboundp;
        int out_of_range_case = 0;

        //        if (lower0 == table.begin() + index1UB) {
        //#ifndef NDEBUG
        //            printf("Warning: val %g too large for row0 lookup with
        //            max of %g. "
        //                   "Setting to a max sbound of %g \n",
        //                   val, *(lower0), para_grid_[colIndex0]);
        //#endif
        //            sbound0 = para_grid_[colIndex0];
        //        } else {
        //            double val0A = *(lower0 - 1);
        //            double val0B = *(lower0);
        //            sbound0 = para_grid_[colIndex0] +
        //                      (val - val0A) / (val0B - val0A) *
        //                      para_spacing_;
        //        }
        // Interpolate first row in the sbound direction
        if (lower1 == table_.begin() + index1UB) {
#ifndef NDEBUG
            printf("Warning: val %g too large for row1 lookup with max of %g. "
                   "Setting sbound max to %g \n",
                   val, *(lower1), para_grid_[colIndexm]);
#endif
            out_of_range_case = 1;
            sboundm = para_grid_[colIndexm]; // Go to end of sbound grid
        } else {
            double valmA = *(lower1 - 1);
            double valmB = *(lower1);
            sboundm = para_grid_[colIndexm] +
                      ((val - valmA) / (valmB - valmA)) * para_spacing_;
        }
        // Interpolate second row in the sbound direction
        if (lower2 == table_.begin() + index2UB) {
#ifndef NDEBUG
            printf("Warning: val %g too large for row2 lookup with max of %g. "
                   "Setting sbound max to %g \n",
                   val, *(lower2), para_grid_[colIndexp]);
#endif
            assert(out_of_range_case != 1); // Something has gone horribly wrong
            out_of_range_case = 2;
            sboundp = para_grid_[colIndexp]; // Go to end of sbound grid
            // return D * sboundm;
        } else {
            double valpA = *(lower2 - 1);
            double valpB = *(lower2);
            sboundp = para_grid_[colIndexp] +
                      ((val - valpA) / (valpB - valpA)) * para_spacing_;
        }
        // printf("sboundm = %f\n", sboundm);
        // printf("sboundp = %f\n", sboundp);

        // Interpolate two sbound values in the dist_perp direction
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
        double UB_val1 = *(table_.begin() + index1UB);
        double UB_val2 = *(table_.begin() + index2UB);
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
        assert(colMax < para_grid_num_);
        // Simple binary search using linear interpolation to find value
        while ((colMax - 1) > colMin) {
            double avg = (colMax + colMin) * .5;
            int colAvg = floor(avg);
            double colVal =
                table_[getTableIndex(rowMin, colAvg)] * (1 - rowFrac)  //
                + table_[getTableIndex(rowMin + 1, colAvg)] * rowFrac; //
            // printf("colAvg = %d, colVal = %f, val = %f\n", colAvg,
            // colVal, val);

            // Accuracy of table is set to 1e-4, use 1e-5 to set colMax
            // value. Otherwise, you will always reach the end of column
            if ((val - colVal) < 1e-5)
                colMax = colAvg;
            else if (colVal < val)
                colMin = colAvg;
        }
        // printf("colMin final = %d\n", colMin);
        // printf("colMax final = %d\n", colMax);
        // Interpolate using function of the line
        double valMin = table_[getTableIndex(rowMin, colMin)] * (1 - rowFrac) //
                        +
                        table_[getTableIndex(rowMin + 1, colMin)] * rowFrac;  //
        double valMax = table_[getTableIndex(rowMin, colMax)] * (1 - rowFrac) //
                        +
                        table_[getTableIndex(rowMin + 1, colMax)] * rowFrac; //

        // Linear interpolation with known values
        double sbound = (val - valMin) / (valMax - valMin) * para_spacing_ +
                        para_grid_[colMin];
        // Make sure sbound is not calculated to be outside range [sMin,
        // sMax]
        if (sbound > sMax)
            return sMax;
        else if (sbound < sMin)
            return sMin;
        else
            return sbound;
    }

    inline int getTableIndex(int row, int col) const {
        assert(row < perp_grid_num_ && row >= 0);
        assert(col < para_grid_num_ && col >= 0);
        return row * para_grid_num_ + col;
    }

  protected:
    /*! \brief Fills the lookup table with values after grid has been
     * created.
     *
     * Table is a flattened 2D array with first index corresponding to
     * different distances head is away from rod and second index is defines
     * the bounds of where that crosslinking head can bind.
     *
     * \return void
     */
    // void FillMatrix() {
    //    table.resize(perp_grid_num_ * para_grid_num_, 0);

    //    //// boost integration parameters
    //    // const int max_depth = 10; // maximum number of interval
    //    splittings
    //    // const double tol = 1e-8;  // maximum relative error

    //    // row major
    //    // i is slow changing, should be perp_grid_
    //    // j is fast changing, should be para_grid_

    //    // IMPORTANT: Table stores dimensionless val
    //    for (int i = 0; i < perp_grid_num_; i++) {
    //        for (int j = 0; j < para_grid_num_; j++) {
    //            const double sbound = para_grid_[j];
    //            assert(!(sbound < 0));
    //            // double result = integral(perp_grid_[i], 0, sbound, M_,
    //            // ell0_);
    //            double result = getIntegralResult(perp_grid_[i], 0,
    //            sbound); table_[i * para_grid_num_ + j] = result;
    //        }
    //    }
    //}

    /**
     * @brief get row index with dist_perp
     *
     * @param dist_perp dimensionless
     * @return int
     */
    inline int getRowIndex(double dist_perp, double &rowFrac) const {
        const double x = (dist_perp - perp_grid_[0]) * perp_spacing_inv_;
        int index = floor(x);
        rowFrac = x - index;
        return index;
    }

    /**
     * @brief Get col index with dist_para
     *
     * @param dist_para dimensionless
     * @return int
     */
    inline int getColIndex(double dist_para, double &colFrac) const {
        const double x = (dist_para - para_grid_[0]) * para_spacing_inv_;
        int index = floor(x);
        colFrac = x - index;
        return index;
    }
};

#endif
