/**
 * @file SylinderNear.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief Essential type for sylinder short range interactions
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef EXAMPLEROD_HPP_
#define EXAMPLEROD_HPP_

//#include "Sylinder.hpp"

//#include "Collision/CollisionCollector.hpp"
//#include "Collision/DCPQuery.hpp"
//#include "FDPS/particle_simulator.hpp"
//#include "Util/EigenDef.hpp"

#include "KMC/kmc.hpp"
#include "KMC/macros.hpp"
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <type_traits>
#include <vector>

#include <mpi.h>
#include <omp.h>

/**
 *  Do not use PS::F64vec3 except in required function interfaces
 *
 */

/**
 * @brief Essential type for sylinder short range interactions
 *
 * Essential Particle Class for FDPS
 */
class ExampleRod {
  public:
    int gid;             ///< global unique id
    double radius;       ///< radius
    double length;       ///< length
    double pos[3];       ///< position
    double direction[3]; ///< direction (unit norm vector)
    double rank;         ///< mpi rank (not necessary but useful)

    /*******************
     *  Get functions  *
     *******************/

    /**
     * @brief Get gid
     *
     * @return int
     */
    int getGid() const { return gid; }

    /**
     * @brief Get center of rod position
     *
     * @return const double*
     */
    const double *getPos() const { return pos; }

    /**
     * @brief Get direction of rod
     *
     * @return const double*
     */
    const double *getDirection() const { return direction; }

    /*******************
     *  Set functions  *
     *******************/

    /**
     * @brief Set center of rod position
     *
     * @return void
     */
    void setPos(const double pos_[]) {
        for (int i = 0; i < 3; ++i) {
            pos[i] = pos_[i];
        }

        return;
    }

    /**
     * @brief Set direction of rod to a unit vector
     *
     * @return
     */
    void setDirection(const double direction_[]) {
        double mag2 = 0;
        for (int i = 0; i < 3; ++i) {
            direction[i] = direction_[i];
            mag2 += direction[i] * direction[i];
        }
        // Make sure direction is indeed a unit vector
        if (ABS(mag2) > 10e-12) {
            double mag = sqrt(mag2);
            for (int i = 0; i < 3; ++i) {
                direction[i] /= mag2;
            }
        }
        return;
    }
};

static_assert(std::is_trivially_copyable<ExampleRod>::value, "");
static_assert(std::is_default_constructible<ExampleRod>::value, "");

#endif
