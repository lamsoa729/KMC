/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : kmc_choose
 * @created     : Saturday Aug 10, 2019 10:45:18 EDT
 */

#ifndef KMC_CHOOSE_HPP

#define KMC_CHOOSE_HPP
//#include "KMC/kmc.hpp"
// class kmc_choose
//{
//    public:
//        kmc_choose ();
//        virtual ~kmc_choose ();
//    private:
//        [> private data <]
//};
//

// template <typename TRod>
inline int choose_kmc_double(double kmc0Prob, double kmc1Prob, double &roll) {
    assert(roll <= 1.0 && roll >= 0.);
    int activated_end = -1;
    // Get probabilities for binding and unbinding (non-normalized)
    double totActivateProb = kmc0Prob + kmc1Prob;
    double passProb = exp(-1. * totActivateProb); // From Poisson process
    // Rescale roll to span unnormalized probabilities
    roll *= (totActivateProb + passProb);
    if (roll < totActivateProb) {
        if (roll < kmc0Prob) {
            activated_end = 0;
            // Renormalize random number to kmc 0 binding
            roll /= kmc0Prob;
        } else {
            activated_end = 1;
            // Shift and renormalize random number to kmc 1 binding
            roll = (roll - kmc0Prob) / kmc1Prob;
        }
    }
#ifndef NDEBUG
    else if (totActivateProb > 0) {
        printf("kmc0 prob: %f, kmc1 prob: %f, no kmc prob: %f,\n", kmc0Prob,
               kmc1Prob, passProb);
    }
#endif
    return activated_end;
}

#endif /* end of include guard KMC_CHOOSE_HPP */

