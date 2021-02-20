#ifndef HELPERS_HPP_
#define HELPERS_HPP_

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "macros.hpp"

/**
 * @brief: Dot product of 3D vectors
 *
 * @param: const double *a
 *       : const double *b
 *
 * @return: double
 */
inline double dot3(const double *a, const double *b) {
    double c = 0;
    for (int i = 0; i < 3; ++i) {
        c += a[i] * b[i];
    }
    return c;
}

/**
 * @brief: Copy one 3D vectors to another
 *
 * @param: const double *a
 *       : const double *b
 *
 * @return: void
 */
inline void copy3(const double *a, double *b) {
    for (int i = 0; i < 3; ++i) {
        b[i] = a[i];
    }
}

/**
 * @brief: Compute the magnitude of a length-3 array.
 *
 * @param: const double *arr3
 *
 * @return: double
 */
inline double magnitude(const double *arr3) {
    double magSqr = SQR(arr3[0]) + SQR(arr3[1]) + SQR(arr3[2]);
    if (magSqr < 0) {
#ifndef NDEBUG
        printf("Warning: magnitude squared (%f) is less than zero. "
               "Setting magnitude to zero.\n", magSqr);
#endif
        return 0;
    }
    return sqrt(magSqr);
}

/**
 * @brief: Find the surface area of the cap of a circle/sphere of radius r1
 * that intersects with a sphere of radius r2 at a center-to-center distance
 * d.
 *
 * @param: const double r1
 * @param: const double r2
 * @param: const double d
 * @param: const int nDim
 *
 * @return: double
 */
/*inline double areaIntersect(double r1, double r2, double d, int nDim) {
    double A = 0;
    if (nDim == 3) {
        2*M_PI*SQR(r1) * (1 + (SQR(r2) - SQR(r1) - SQR(d))/(2*r1*d));
        if (A > 4*M_PI*SQR(r1)) {
#ifndef NDEBUG
        printf("Warning: intersect area %f greater than possible", A);
#endif
*/


/**
 * @brief: Find the minimum distance between a point and a line segment
 *
 * @param: const double *point
 * @param: const double *minus
 * @param: const double *plus
 * @param: double *pointPerp
 *
 * @return: double
 */
inline double dist_point_seg(const double *point, const double *minus,
                             const double *plus, double *pointPerp) {
    // The direction vector is not unit length.  The normalization is deferred
    // until it is needed.
    double direction[3], diff[3];
    for (int i = 0; i < 3; ++i) {
        direction[i] = plus[i] - minus[i];
        diff[i] = point[i] - plus[i];
    }
    double t = dot3(direction, diff);
    if (t >= 0) {
        copy3(plus, pointPerp);
    } else {
        for (int i = 0; i < 3; ++i) {
            diff[i] = point[i] - minus[i];
        }
        t = dot3(direction, diff);
        if (t <= 0) {
            copy3(minus, pointPerp);
        } else {
            double sqrLength = dot3(direction, direction);
            if (sqrLength > 0) {
                t /= sqrLength;
                for (int i = 0; i < 3; ++i) {
                    pointPerp[i] = minus[i] + t * direction[i];
                }
            } else {
                copy3(minus, pointPerp);
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        diff[i] = point[i] - pointPerp[i];
    }
    return sqrt(dot3(diff, diff));
};

template <typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif /* HELPERS_HPP */
