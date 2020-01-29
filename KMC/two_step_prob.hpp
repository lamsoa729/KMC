/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : two_step_prob
 * @created     : Wednesday Jan 29, 2020 14:05:35 MST
 */

#ifndef TWO_STEP_PROB_HPP

#define TWO_STEP_PROB_HPP
#include <boost/math/tools/roots.hpp>
#include <cmath>

// Solve the for the time that maximizes the probability of two events occuring
template <class T>
class MaxTimeOfTwoStepProbFunctor {
    // Constructor just stores values find root of function.
    MaxTimeOfTwoStepProbFunctor(T const &rate1, T const &rate2,
                                T const &total_time, T const &tolerance)
        : k_ab(rate1), k_bc(rate2), tot_time(total_time), eps(tolerance) {}
    T operator()(T const &t) {
        T fx = k_ab * (exp(k_bc * (tot_time - t)) - k_bc * exp(k_ab * t)) - eps;
        return fx;
    }

  private:
    T k_ab;
    T k_bc;
    T tot_time;
    T eps;
};

template <class T>
T MaxTimeOfTwoStepProb(T const &rate1, T const &rate2, T const &total_time,
                       T const &tolerance) {

    T guess = .5 * total_time;
    T factor = 2; // How big steps to take when searching.

    const boost::uintmax_t maxit = 20; // Limit to maximum iterations.
    boost::uintmax_t it =
        maxit; // Initally our chosen max iterations, but updated with actual.
    bool is_rising =
        true; // So if result if guess^3 is too low, then try increasing guess.
    int digits = std::numeric_limits<T>::digits; // Maximum possible binary
                                                 // digits accuracy for type T.
    // Some fraction of digits is used to control how accurate o try to make the
    // result.
    int get_digits =
        digits - 3; // We have to have a non-zero interval at each step, so
                    // maximum accuracy is digits - 1.  But we also have to
                    // allow for inaccuracy in f(x), otherwise the last few
                    // iterations just thrash around.
    boost::math::tools::eps_tolerance<T> rel_tol(
        get_digits); // Set the tolerance.
    std::pair<T, T> t = boost::math::tools::bracket_and_solve_root(
        MaxTimeOfTwoStepProbFunctor<T>(rate1, rate2, total_time, tolerance),
        guess, factor, is_rising, rel_tol, it);
    return t.first + (t.second - t.first) * .5;
};

#endif /* end of include guard TWO_STEP_PROB_HPP */

