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

#include "integrals.hpp"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

/**
 * @brief 2D Lookup Table for the dimensionless functions
 *
 * Init() is not thread-safe
 * Lookup() and ReverseLookup() are thread safe
 *
 * Both interface functions and internal caculations are dimensionless
 */
class LookupTable {
  public:
    double D; ///< dimensional tubule Diameter, length scale

    double ell0; ///< \ell_0/D, dimensionless
    double M;    ///< (1-\lambda)\kappa\beta/2 * D^2, dimensionless

    int distPerpGridNumber;
    double distPerpGridSpacingInv;    ///< dimensionless
    std::vector<double> distPerpGrid; ///< dimensionless vertical direction

    int sboundGridNumber;
    double sboundGridSpacingInv;    ///< dimensionless
    std::vector<double> sboundGrid; ///< dimensionless horizontal direciton

    std::vector<double> table; ///< the 2D matrix of dimensionless values

  public:
    LookupTable() = default;
    ~LookupTable() = default;

    /**
     * @brief Get the TubuleDiameter length scale
     *
     * @return double
     */
    double getTubuleDiameter() const { return D; }

    /**
     * @brief Get the dimensionless sbound for tabulation
     *
     * @return double
     */
    double getNonDsbound() const { return sboundGrid.back(); }

    /**
     * @brief Initialize the lookup table with alpha
     *
     * @param M the factor of exponential (1-\lambda)\kappa\beta/2
     * @param freeLength the protein free length
     * @param tubuleD the MicroTubule diameter
     */
    void Init(double M_, double freeLength, double tubuleD);

    /*********************
     *  output LookupTable
     *********************/

    /**
     * @brief Print Table to std::cout
     * Currently 2D only
     */
    void PrintTable() const;

    /*********************
     *  Look Up
     *********************/

    /** TODO: test this
     * @brief Lookup two points
     *
     * @param distPerp dimensionless
     * @param sbound dimensionless
     * @return dimensionless integration value
     */
    double Lookup(double distPerp, double sbound) const;

    /******************
     * Invert Lookup
     ******************/

    /** TODO: test this
     * @brief Reverse lookup
     *
     * @param distPerp dimensionless
     * @param val dimensionless
     * @return sbound dimensionless
     */
    double ReverseLookup(double distPerp, const double val) const;

  private:
    void FillMatrix();

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

    inline int getTableIndex(int row, int col) const {
        assert(row < distPerpGridNumber && row >= 0);
        assert(col < sboundGridNumber && col >= 0);
        return row * sboundGridNumber + col;
    }
};

#endif
