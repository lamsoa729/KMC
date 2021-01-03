/**
 * @file ExampleSphere.hpp
 * @author Modified from code by wenyan4work@gmail.com
 * @brief Generic sphere type
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef EXAMPLESPHERE_HPP_
#define EXAMPLESPHERE_HPP_

#include "KMC/kmc.hpp"
#include "KMC/macros.hpp"
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <type_traits>
#include <vector>

template <class TSphere>
TSphere MockSphere(int id) {
    TSphere sphere;
    sphere.gid = id;
    sphere.rank = 0;
    sphere.radius = .5;
    for (int i = 0; i < 3; ++i) sphere.pos[i] = 0;
    return sphere;
}

/**
 * @brief Type definition for generic sphere
 *
 * 
 */
class ExampleSphere {
  public:
    int gid;             ///< global unique id
    double radius;       ///< radius
    double pos[3];       ///< position
    double rank;         ///< mpi rank (not implemented now)

    /*******************
     *  Get functions  *
     *******************/
    int getGid() const { return gid; }
    const double *getPos() const { return pos; }
    void setPos(const double pos_[]) {
        for (int i = 0; i < 3; ++i) {
            pos[i] = pos_[i];
        }
        return;
    }
};

static_assert(std::is_trivially_copyable<ExampleSphere>::value, "");
static_assert(std::is_default_constructible<ExampleSphere>::value, "");

#endif
