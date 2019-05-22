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

// dimensionless integral gaussian integral
double integral(double lm, double sbound0, double sbound1, double M,
                double ell0);

#endif /* INTEGRALS_HPP_ */
