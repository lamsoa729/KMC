/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : lut_filler_base
 * @created     : Friday Feb 07, 2020 15:46:13 MST
 */

#ifndef LUT_FILLER_EDEP_HPP

#define LUT_FILLER_EDEP_HPP
#include "integrals.hpp"
#include "lut_filler_base.hpp"

#include <cassert>
#include <vector>

class LUTFillerEdep : public LUTFillerBase {
  protected:
    double exp_fact_;
    double rest_length_;
    double D_;

  public:
    static constexpr double small_ = 1e-4;

    LUTFillerEdep(double dist_para_grid_num, double dist_perp_grid_num)
        : LUTFillerBase(dist_para_grid_num, dist_perp_grid_num) {}

    // virtual ~LUTFillerEdep();

    void Init(double exp_fact, double rest_length, double D);
    // void Init(double exp_fact, double rest_length, double D) {
    //    D_ = D;
    //    exp_fact_ = exp_fact * D_ * D_;
    //    rest_length_ = rest_length / D_;
    //    LUTFillerBase::Init();
    //}

    double getUpperBound() const;
    // double getUpperBound() const {
    //    return sqrt(-log(small_) / exp_fact_) + rest_length_;
    //}
    double getIntegralResult(double dist_perp, double dist_para_l,
                             double dist_para_u) const;
    // double getIntegralResult(double dist_perp, double dist_para_l,
    //                         double dist_para_u) const {
    //    return integral(dist_perp, dist_para_l, dist_para_u, exp_fact_,
    //                    rest_length_);
    //}
    double getBindingVolume() const;
    // double getBindingVolume() const {
    //    return bind_vol_integral(upper_bound_, exp_fact_, rest_length_);
    //}
    inline const double getDiameter() const { return D_; }
};

void LUTFillerEdep::Init(double exp_fact, double rest_length, double D) {
    D_ = D;
    exp_fact_ = exp_fact * D_ * D_;
    rest_length_ = rest_length / D_;
    LUTFillerBase::Init();
}
double LUTFillerEdep::getUpperBound() const {
    return sqrt(-log(small_) / exp_fact_) + rest_length_;
}
double LUTFillerEdep::getIntegralResult(double dist_perp, double dist_para_l,
                                        double dist_para_u) const {
    return integral(dist_perp, dist_para_l, dist_para_u, exp_fact_,
                    rest_length_);
}
double LUTFillerEdep::getBindingVolume() const {
    return bind_vol_integral(upper_bound_, exp_fact_, rest_length_);
}

#endif /* end of include guard LUT_FILLER_EDEP_HPP */
