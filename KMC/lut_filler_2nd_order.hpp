
/**
 * @file lut_filler_2nd_order.hpp
 * @author Adam Lamson
 * @brief Class to fill lookup table with
 * @version 0.1
 * @date 2020-02-10
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef LUT_FILLER_2ND_ORDER_HPP_
#define LUT_FILLER_2ND_ORDER_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "integrals.hpp"
#include "lut_filler_edep.hpp"
#include "macros.hpp"
/**
 * @brief 2D Lookup Table for the dimensionless functions
 *
 *
 * Both interface functions and internal caculations are dimensionless
 * 
 * Only differ between this and lut_filler_edep is the integral to calculate 
 * the lookup table adds a slight correction which is more accurate when filaments
 * are mostly parallel.
 * 
 */
class LUTFiller2ndOrder : public LUTFillerEdep {

  public:
    static constexpr double small_ = 1e-4;
    LUTFiller2ndOrder(double dist_para_grid_num, double dist_perp_grid_num)
        : LUTFillerEdep(dist_para_grid_num, dist_perp_grid_num) {}

    /*! \brief Initialize lookup table filler for asymmetric spring
     *
     * Stores precalculated values of the integral
     * To do this, grid spacing are calculated based on the value of exp_fact_
     * and freeLength. The upper bound lUB is calculated so that the maximum
     * change at the end of the integral is less than 10^{-4}. Everything is
     * non-dimensionalized by rodD.
     *
     * \param exp_fact Exponential prefactor when spring is stretched with the
     * dimensions of L^{-2}.
     * \param freeLength Rest length of binding crosslink.
     * \param rodD Diameter of rod crosslink is binding to.
     *
     * \return void
     */
    void Init(double exp_fact, double rest_length, double D) {
        LUTFillerEdep::Init(exp_fact, rest_length, D);
    }

    double getIntegralResult(double dist_perp, double dist_para_l,
                             double dist_para_u) const {
        return edep_second_order_integral(dist_perp, dist_para_l, dist_para_u, exp_fact_,
                        rest_length_);
    }
    double getUpperBound() const {
        return sqrt(-log(small_) / exp_fact_) + rest_length_ + 5.;
    }
};

#endif
