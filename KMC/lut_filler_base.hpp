/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : lut_filler_base
 * @created     : Friday Feb 07, 2020 15:46:13 MST
 */

#ifndef LUT_FILLER_BASE_HPP

#define LUT_FILLER_BASE_HPP
#include "integrals.hpp"
#include <array>
#include <cassert>

template <typename T>
class LUTFillerBase {
  protected:
    const int dist_para_grid_num_;
    const int dist_perp_grid_num_;

    std::vector<double> dist_para_grid_;
    std::vector<double> dist_perp_grid_;

    double D_;

    double upper_bound_ = -1;

  public:
    LUTFillerBase(double dist_para_grid_num, double dist_perp_grid_num)
        : dist_para_grid_num_(dist_para_grid_num),
          dist_perp_grid_num_(dist_perp_grid_num) {}

    virtual ~LUTFillerBase();

    void Init() {
        upper_bound_ = getUpperBound();
        fillDistParaGrid(dist_para_grid_);
        fillDistPerpGrid(dist_perp_grid_);
    }

    virtual double getUpperBound() const;
    virtual double getBindingVolume() const;
    virtual double getIntegralResult(double dist_perp, double dist_para_l,
                                     double dist_para_u) const;

    inline const double getDistParaGridSpacing() const {
        assert(dist_para_grid_num_ > 1);
        return (upper_bound_) / (dist_para_grid_num_ - 1);
    }
    inline const double getDistPerpGridSpacing() const {
        assert(dist_perp_grid_num_ > 1);
        return (upper_bound_) / (dist_perp_grid_num_ - 1);
    }
    inline const double getDiameter() const { return D_; }

    void fillDistParaGrid(std::vector<double> &dist_para_grid) {
        double spacing = getDistParaGridSpacing();
        dist_para_grid.resize(dist_para_grid_num_);
        for (int i = 0; i < dist_para_grid_num_; i++) {
            dist_para_grid[i] = i * spacing;
        }
    }
    void fillDistPerpGrid(std::vector<double> &dist_perp_grid) {
        double spacing = getDistPerpGridSpacing();
        dist_perp_grid.resize(dist_perp_grid_num_);
        for (int i = 0; i < dist_perp_grid_num_; i++) {
            dist_perp_grid[i] = i * spacing;
        }
    }

    void fillMatrix(std::vector<double> &table) const {
        table.resize(dist_para_grid_num_ * dist_perp_grid_num_, 0);

        //// boost integration parameters
        // const int max_depth = 10; // maximum number of interval splittings
        // const double tol = 1e-8;  // maximum relative error

        // row major
        // i is slow changing, should be distPerpGrid
        // j is fast changing, should be sboundGrid

        // IMPORTANT: Table stores dimensionless val
        for (int i = 0; i < dist_perp_grid_num_; i++) {
            for (int j = 0; j < dist_para_grid_num_; j++) {
                const double s = dist_para_grid_[j];
                assert(!(s < 0));
                double result = getIntegralResult(dist_perp_grid_[i], 0,
                                                  dist_para_grid_[j]);
                table[i * dist_para_grid_num_ + j] = result;
            }
        }
    }
};

#endif /* end of include guard LUT_FILLER_BASE_HPP */
