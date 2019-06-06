#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_
/**
 * @file integrals.hpp
 * @author Adam Lamson
 * @brief Probility Density Function integrals for KMC algorithm. Used to create
 * lookup tables
 * @version 0.1
 * @date 2019-04-15
 *
 * @copyright Copyright (c) 2019
 *
 */

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>

/*! \brief Integrate exponential factor with the form of
 * e^{-M * (\sqrt{ s^2 + lm^2}-ell0 -1)^2}
 * from sbound0 to sbound1  with respect to the variable s.
 *
 * \param lm Physically, this is the perpendicular distance above rod
 * \param sbound lowerr limit of integral
 * \param sbound Upper limit of integral
 * \param M exponential constant factor. Physically, this is the product of
 (1-load_sensitivity)*spring_const/(k_B * Temperature)
 * \param ell0 Shift of the integrands mean. Physically, protein rest length
 * \return result The value of the integration

 */
inline double integral(double lm, double sbound0, double sbound1, double M,
                       double ell0) {
    if (sbound0 >= sbound1) {
        return 0;
    }
    auto integrand = [&](double s) {
        // lambda capture variabls ell0 and M
        const double exponent = sqrt(s * s + lm * lm) - ell0 - 1;
        return exp(-M * exponent * exponent);
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

#endif /* INTEGRALS_HPP_ */
