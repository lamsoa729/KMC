/**
 * @file fdep_lookup_table.hpp
 * @author Wen Yan
 * @brief 2D lookup table for probability CDF calculations
 * @version 0.1
 * @date 2019-03-08
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef FDEP_LOOKUP_TABLE_HPP_
#define FDEP_LOOKUP_TABLE_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "integrals.hpp"
#include "lookup_table.hpp"
#include "macros.hpp"
/**
 * @brief 2D Lookup Table for the dimensionless functions
 *
 * Init() is not thread-safe
 * Lookup() and ReverseLookup() are thread safe
 *
 * Both interface functions and internal caculations are dimensionless
 */
class FdepLookupTable : public LookupTable {
  private:
    double e_fact_;      ///< \lambda, dimensionless
    double fdep_length_; ///< xc/ D, dimensionless

  public:
    FdepLookupTable() = default;
    ~FdepLookupTable() = default;

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
    void Init(double M, double e_fact, double fdep_length, double freeLength,
              double rodD) {
        // constexpr double small = 1e-4;
        // setup length scale and dimensionless lengths
        D_ = rodD;
        ell0_ = freeLength / rodD;
        fdep_length_ = fdep_length / rodD;
        e_fact_ = e_fact;
        M_ = M * rodD * rodD;
        assert(D_ > 0);
        assert(0 <= e_fact_ && e_fact <= 1);
        assert(M_ > 0);

        // Truncate the integration when integrand < SMALL
        // interation table in dimensionless lengths
        // Solve ln(small)/M = .5(1-e_fact)(r-ell0)^2 -
        // 2*fdep_length*(r-ell0) for r
        lUB_ = ell0_ + ((sqrt(SQR(fdep_length_) -
                              (2. * (1. - e_fact_) * log(small_) / M_)) +
                         fdep_length_) /
                        (1. - e_fact_));
        assert(lUB_ > 0);

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
    void calcBindVol() {
        bind_vol_ =
            fdep_bind_vol_integral(lUB_, M_, e_fact_, fdep_length_, ell0_);
        printf("bind_vol_ = %f\n", bind_vol_);
    }
};

#endif
