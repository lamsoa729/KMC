/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : lut_filler_base
 * @created     : Friday Feb 07, 2020 15:46:13 MST
 */

#ifndef LUT_FILLER_EDEP_HPP

#define LUT_FILLER_EDEP_HPP
#include "integrals.hpp"
#include "lut_filler.hpp"

#include <cassert>
#include <vector>

class LUTFillerEdep : public LUTFiller {
  protected:
    double exp_fact_;
    double rest_length_;

  public:
    static constexpr double small_ = 1e-4;

    LUTFillerEdep(double dist_para_grid_num, double dist_perp_grid_num)
        : LUTFiller(dist_para_grid_num, dist_perp_grid_num) {}

    // virtual ~LUTFillerEdep();

    /*! \brief Initialize lookup table filler
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
    void Init(double exp_fact, double rest_length, double D) {
        length_scale_ = D;
        exp_fact_ = exp_fact * length_scale_ * length_scale_;
        rest_length_ = rest_length / length_scale_;
        LUTFiller::Init();
    }

    // Getter functions for Boltzmann variables
    double getExpFact1() const { return exp_fact_; }
    double getExpFact2() const { return exp_fact_; }
    double getEFact() const { return 0; } // in exp_fact
    double getFDepLength() const { return 0; }
    double getRestLength() const { return rest_length_; }

    // Calculates the Boltzmann factor for point-like object
    inline double calcBoltzmann(double dist_cent) const {
        /* TODO: Make dimensionful unit test for this <30-03-21, ARL> */
        double r = dist_cent / length_scale_;
        return exp(-exp_fact_ * SQR(r - rest_length_));
    }

    // double getUpperBound() const;
    double getUpperBound() const {
        return sqrt(-log(small_) / exp_fact_) + rest_length_;
    }
    // double getIntegralResult(double dist_perp, double dist_para_l,
    //                         double dist_para_u) const;
    double getIntegralResult(double dist_perp, double dist_para_l,
                             double dist_para_u) const {
        return integral(dist_perp, dist_para_l, dist_para_u, exp_fact_,
                        rest_length_);
    }
    // double getBindingVolume() const;
    double getBindingVolume() const {
        return bind_vol_integral(upper_bound_, exp_fact_, rest_length_);
    }
};

#endif /* end of include guard LUT_FILLER_EDEP_HPP */
