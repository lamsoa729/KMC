/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : two_step_prob
 * @created     : Wednesday Jan 29, 2020 14:05:35 MST
 */

#ifndef TWO_STEP_PROB_HPP

#define TWO_STEP_PROB_HPP
#include <boost/math/tools/roots.hpp>
#include <cassert>
#include <cmath>

// Solve the for the time that maximizes the probability of two events occuring
template <class T>
class MaxTimeOfTwoStepProbFunctor {
  public:
    // Constructor just stores values find root of function.
    MaxTimeOfTwoStepProbFunctor(T const &rate1, T const &rate2,
                                T const &total_time)
        : k_ab(rate1), k_bc(rate2), tot_time(total_time) {}

    T operator()(T const &t) {
        T fx = k_ab * (exp(k_bc * (tot_time - t)) - 1.) -
               k_bc * (exp(k_ab * t) - 1.);
        return fx;
    }

  private:
    T k_ab;
    T k_bc;
    T tot_time;
};

template <class T>
T max_time_of_two_step_prob(T const &rate1, T const &rate2,
                            T const &total_time) {

    T guess = .5 * total_time;
    T factor = 2; // How big steps to take when searching.

    const boost::uintmax_t maxit = 20; // Limit to maximum iterations.
    boost::uintmax_t it =
        maxit; // Initally our chosen max iterations, but updated with actual.
    bool is_rising = false; // If guess is too low, try increasing guess.
    int digits = std::numeric_limits<T>::digits; // Maximum possible binary
                                                 // digits accuracy for type T.
    // Fraction of digits used to control accuracy of result.
    int get_digits =
        digits - 3; // We have to have a non-zero interval at each step, so
                    // maximum accuracy is digits - 1.  But we also have to
                    // allow for inaccuracy in f(x), otherwise the last few
                    // iterations just thrash around.
    boost::math::tools::eps_tolerance<T> rel_tol( get_digits); // Set tolerance.
    std::pair<T, T> t = boost::math::tools::bracket_and_solve_root(
        MaxTimeOfTwoStepProbFunctor<T>(rate1, rate2, total_time), guess, factor,
        is_rising, rel_tol, it);
    return t.first + (t.second - t.first) * .5; // Zero occurs between values
};

template <class T>
T two_step_prob(T const &rate1, T const &rate2, T const &total_time,
                T const &t) {
    assert(total_time >= t);
    T prob = (1. - exp(-1. * rate1 * t)) *
             (1. - exp(-1. * rate2 * (total_time - t)));
    return prob;
};

template <class T>
T two_step_max_prob(T const &rate1, T const &rate2, T const &total_time) {
    if (rate1 == 0. && rate2 == 0.)
        return 0.;
    if (rate1 == 0.)
        return 1. - exp(-1. * rate2 * total_time);
    if (rate2 == 0.)
        return 1. - exp(-1. * rate1 * total_time);

    const T t_max = max_time_of_two_step_prob(rate1, rate2, total_time);
    if (t_max >= total_time || std::isnan(t_max)){
        printf("WARNING!!!: Time of max probability of two bind steps occuring " 
               "is [%f] whereas a timestep is [%f]. " 
               "Setting probability using largest rate.\n", t_max, total_time); 
        const T rate = rate1 > rate2 ? rate1 : rate2;
        return 1. - exp(-1. * rate * total_time);
    }
    return two_step_prob(rate1, rate2, total_time, t_max);
};

#endif /* end of include guard TWO_STEP_PROB_HPP */

