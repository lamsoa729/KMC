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
 * e^{-M * (\sqrt{ s^2 + lm^2} - ell0)^2}
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
        const double exponent = sqrt(s * s + lm * lm) - ell0;
        return exp(-M * exponent * exponent);
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

inline double bind_vol_integral(double sbound, double M, double ell0) {
    assert(sbound > 0);
    auto integrand = [&](double s) {
        // lambda capture variabls ell0 and M
        const double exponent = s - ell0;
        return s * s * exp(-M * exponent * exponent);
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, 0, sbound, 10, 1e-6, &error);
    return 4. * M_PI * result;
}

/*! \brief Integrate exponential factor with the form of
 * e^{-M * [ (1-e_fact)\sqrt{s^2 + lm^2} - ell0)^2 -
 *          fdep_length * (\sqrt{ s^2 + lm^2} - ell0) ] }
 * from sbound0 to sbound1  with respect to the variable s.
 *
 * \param lm Physically, this is the perpendicular distance above rod
 * \param sbound lowerr limit of integral
 * \param sbound Upper limit of integral
 * \param M exponential constant factor. Physically, this is the product of
 spring_const/(k_B * Temperature)
 * \param e_fact energy(load) sensitivity to unbinding.
 * \param fdep_length Characteristic length for force dependent unbinding.
 * \param ell0 Shift of the integrands mean. Physically, protein rest length
 * \return result The value of the integration

 */
inline double fdep_integral(double lm, double sbound0, double sbound1, double M,
                            double e_fact, double fdep_length, double ell0) {
    if (sbound0 >= sbound1) {
        return 0;
    }
    auto integrand = [&](double s) {
        const double rprime = sqrt(s * s + lm * lm) - ell0;
        const double energy_term = .5 * (1. - e_fact) * rprime * rprime;
        const double force_term = fdep_length * rprime;
        return exp(-M * (energy_term - force_term));
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

inline double fdep_bind_vol_integral(double sbound, double M, double e_fact,
                                     double fdep_length, double ell0) {
    assert(sbound > 0);
    auto integrand = [&](double s) {
        const double rprime = s - ell0;
        const double energy_term = .5 * (1. - e_fact) * rprime * rprime;
        const double force_term = fdep_length * rprime;
        return exp(-M * (energy_term - force_term));
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, 0, sbound, 10, 1e-6, &error);
    return 4. * M_PI * result;
}

/*! \brief Integrate exponential factor with the form of
 * e^{-M * [ (1-e_fact)\sqrt{s^2 + lm^2} - ell0)^2 -
 *          fdep_length * (\sqrt{ s^2 + lm^2} - ell0) ] }
 * from sbound0 to sbound1  with respect to the variable s.
 *
 * \param lm Physically, this is the perpendicular distance above rod
 * \param sbound lower limit of integral
 * \param sbound Upper limit of integral
 * \param M1 exponential constant factor when spring is compressed.
 *    Physically, this is the product of spring_const_1/(k_B * Temperature)
 * \param M2 exponential constant factor when spring is stretched.
 *    Physically, this is the product of spring_const_2/(k_B * Temperature)
 * \param e_fact energy(load) sensitivity to unbinding.
 * \param fdep_length Characteristic length for force dependent unbinding.
 * \param ell0 Shift of the integrands mean. Physically, protein rest length
 * \return result The value of the integration

 */
inline double asym_integral(double lm, double sbound0, double sbound1,
                            double M1, double M2, double e_fact,
                            double fdep_length, double ell0) {
    if (sbound0 >= sbound1) {
        return 0;
    }
    auto integrand = [&](double s) {
        const double rprime = sqrt(s * s + lm * lm) - ell0;
        const double energy_term = .5 * (1. - e_fact) * rprime * rprime;
        const double force_term = fdep_length * rprime;
        const double M = rprime < 0. ? M1 : M2;
        return exp(-M * (energy_term - force_term));
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

inline double asym_bind_vol_integral(double sbound, double M1, double M2,
                                     double e_fact, double fdep_length,
                                     double ell0) {
    assert(sbound > 0);
    auto integrand = [&](double s) {
        const double rprime = s - ell0;
        const double energy_term = .5 * (1. - e_fact) * rprime * rprime;
        const double force_term = fdep_length * rprime;
        const double M = rprime < 0. ? M1 : M2;
        return exp(-M * (energy_term - force_term));
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, 0, sbound, 10, 1e-6, &error);
    return 4. * M_PI * result;
}

/*! \brief Integrate exponential factor with the form of
 * e^{-M * ( (lm-1)\sqrt{ 1 +(s/lm)^2} - ell0)^2}
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
inline double edep_second_order_integral(double lm, double sbound0, double sbound1, double M,
                       double ell0) {
    if (sbound0 >= sbound1) {
        return 0;
    }
    //TODO add exception for when lm < 1 ?
    auto integrand = [&](double s) {
        const double exponent = ((lm - 1.) * sqrt(1. + (s*s/(lm*lm)) )) - ell0;
        return exp(-M * exponent * exponent);
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

// inline double bind_vol_integral(double sbound, double M, double ell0) {
//     assert(sbound > 0);
//     auto integrand = [&](double s) {
//         // lambda capture variabls ell0 and M
//         const double exponent = s - ell0;
//         return s * s * exp(-M * exponent * exponent);
//     };
//     double error = 0;
//     double result =
//         boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
//             integrand, 0, sbound, 10, 1e-6, &error);
//     return 4. * M_PI * result;
// }

#endif /* INTEGRALS_HPP_ */
