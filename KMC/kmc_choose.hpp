/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : kmc_choose
 * @created     : Saturday Aug 10, 2019 10:45:18 EDT
 */

#ifndef KMC_CHOOSE_HPP

#define KMC_CHOOSE_HPP
#include <cassert>
#include <math.h>
#include <stdio.h>

inline int rescale_roll(const double kmc0Prob, const double kmc1Prob,
                        double &roll) {
    int activated_end = -1;
    if (roll < kmc0Prob) {
        activated_end = 0;
        // Renormalize random number to kmc 0 binding
        roll /= kmc0Prob;
    } else {
        activated_end = 1;
        // Shift and renormalize random number to kmc 1 binding
        roll = (roll - kmc0Prob) / kmc1Prob;
    }
    return activated_end;
}

inline int choose_kmc_double(double kmc0Prob, double kmc1Prob, double &roll) {
    assert(roll <= 1.0 && roll >= 0.);
    assert(kmc0Prob >= 0);
    assert(kmc1Prob >= 0);
    int activated_end = -1;
    if (kmc0Prob == HUGE_VAL && kmc1Prob == HUGE_VAL) {
        activated_end = rescale_roll(.5, .5, roll);
    } else if (kmc0Prob == HUGE_VAL) {
        activated_end = 0;
    } else if (kmc1Prob == HUGE_VAL) {
        activated_end = 1;
    } else {
        // Get probabilities for binding and unbinding (non-normalized)
        double totActivateProb = kmc0Prob + kmc1Prob;
        double passProb = exp(-1. * totActivateProb); // From Poisson process
        // Rescale probabilities to add to 1 so that roll samples all events
        kmc0Prob /= (totActivateProb + passProb);
        kmc1Prob /= (totActivateProb + passProb);
        totActivateProb /= (totActivateProb + passProb);

#ifndef NDEBUG
        if (totActivateProb > .3)
            printf("totActivateProb: %f, passProb: %f \n", totActivateProb,
                   passProb);
#endif
        if (roll < totActivateProb) {
            activated_end = rescale_roll(kmc0Prob, kmc1Prob, roll);
#ifndef NDEBUG
            printf("kmc0 prob: %f, kmc1 prob: %f, no kmc prob: %f \n", kmc0Prob,
                   kmc1Prob, passProb);
#endif
        }
    }
    return activated_end;
}

#endif /* end of include guard KMC_CHOOSE_HPP */

