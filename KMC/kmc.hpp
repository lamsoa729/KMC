#ifndef KMC_HPP
#define KMC_HPP

#include "KMC/helpers.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/macros.hpp"

#include <array>
#include <cassert>
#include <cmath>

template <typename TRod>
class KMC {
  private:
    // System parameters
    const double dt_; // Time step

    // Probabilities
    double prob_tot_ = 0;
    double r_cutoff_;

    std::vector<double> rods_probs_; ///< binding probabilities

    // Spatial variables
    double pos_[3];                   ///< position of motor head
    std::vector<double> distMinArr_;  ///< min dist to rod segment
    std::vector<double> distPerpArr_; ///< min (perp) dist to rod line
    std::vector<double> muArr_; ///< Closest points of proteins to rods along
                                ///< the rods axis with  respect to rod center
                                ///< and direction.
    // std::vector<double> &bindFactorsArr_; ///< Head binding factors

    std::vector<std::pair<double, double>> lims_;
    ///< Bounds of binding positions on rods with respect to closest
    ///< point of rod axis. Used in lookup table value acquisition.

    const LookupTable *LUTablePtr_; ///< lookup table

  public:
    /******************
     *  Constructors  *
     ******************/
    // Constructor for unbinding ends
    KMC(const double *pos, const double dt) : r_cutoff_(-1), dt_(dt) {
        setPos(pos);
    }

    // Constructor for binding ends with diffusion modeled
    // TODO: create unit test for this
    KMC(const double *pos, const int Npj, const double r_cutoff,
        const double DiffConst, const double dt)
        : dt_(dt), LUTablePtr_(nullptr) {
        setPos(pos);

        // Find average diffusion distance and compare to given r_cutoff
        double avg_dist = getDiffRadius(diffConst);
        r_cutoff_ = avg_dist > r_cutoff ? avg_dist : r_cutoff;

        rods_probs_.resize(Npj, 0);
        distMinArr_.resize(Npj, 0);
        distPerpArr_.resize(Npj, 0);
        muArr_.resize(Npj, 0);
        lims_.resize(Npj);
    }

    // Constructor for S->D end binding without lookup tables
    KMC(const double *pos, const int Npj, const double r_cutoff,
        const double dt)
        : r_cutoff_(r_cutoff), dt_(dt), LUTablePtr_(nullptr) {
        setPos(pos);
        rods_probs_.resize(Npj, 0);
        distMinArr_.resize(Npj, 0);
        distPerpArr_.resize(Npj, 0);
        muArr_.resize(Npj, 0);
        lims_.resize(Npj);
    }

    // Constructor for binding ends with lookup tables
    KMC(const double *pos, const int Npj, const double r_cutoff,
        const double dt, const LookupTable *LUTablePtr)
        : r_cutoff_(r_cutoff), dt_(dt), LUTablePtr_(LUTablePtr) {
        setPos(pos);
        rods_probs_.resize(Npj, 0);
        distMinArr_.resize(Npj, 0);
        distPerpArr_.resize(Npj, 0);
        muArr_.resize(Npj, 0);
        lims_.resize(Npj);
    }

    /*************************************
     *  Calculate probability functions  *
     *************************************/

    void CalcTotProbsUS(const TRod *const *rods,
                        const std::vector<int> &uniqueFlagJ,
                        const std::vector<double> &bindFactors);

    double CalcProbUS(const int j_bond, const TRod &rod,
                      const double bindFactor);

    void CalcProbSU(const double unbindFactor);

    void CalcTotProbsSD(const TRod *const *rods,
                        const std::vector<int> &uniqueFlagJ, const int boundID,
                        const double lambda, const double kappa,
                        const double beta, const double restLen,
                        const std::vector<double> &bindFactors);

    double LUCalcProbSD(const int j_rod, const TRod &rod,
                        const double bindFactor);

    double CalcProbSD(const int j_rod, const TRod &rod, const double lambda,
                      const double kappa, const double beta,
                      const double restLen, const double bindFactor);

    void CalcProbDS(const double unbindFactor);

    void UpdateRodDistArr(const int j_bond, const TRod &rod);

    double getTotProb() { return prob_tot_; }

    int whichRodBindUS(const TRod *const *rods, double &bindPos, double roll);

    void whereUnbindSU(double R, double diffConst, double rollVec[3],
                       double pos[3]);

    int whichRodBindSD(double &bindPos, double roll);

    double RandomBindPosSD(int j_bond, double roll);

    void whereUnbindDS(const double boundPos[], double pos[]);

    /*******************
     *  Get functions  *
     *******************/

    double getDiffRadius(const double diffConst) {
        return sqrt(6.0 * diffConst * dt_);
    }

    double getMu(const int j_bond) { return muArr_[j_bond]; }

    double getDistMin(const int j_bond) { return distMinArr_[j_bond]; }

    double getDistPerp(const int j_bond) { return distPerpArr_[j_bond]; }

    double getProbs(const int j_bond) { return rods_probs_[j_bond]; }

    void setPos(const double *pos) {
        for (int i = 0; i < 3; ++i) {
            pos_[i] = pos[i];
        }
    }

    virtual ~KMC() {}
};

/*! \brief Function to update the minimum distance of KMC object to rod
 * (distMinArr_), perpendicular distance from rod axis (distPerpArr_), as well
 * as the point located on the rod from rods center where
 * this minimum distance occurs (muArr_).
 *
 *
 * \param j_rod KMC index of rod used in muArr_ and distPerpArr_
 * \param &rod Rod data structure containing center location, length,
 *        and unit direction vector.
 * \return void, Updates muArr_[j_rod], distMinArr_[j_rod], and
 * distPerpArr_[j_rod] in KMC object.
 */
template <typename TRod>
void KMC<TRod>::UpdateRodDistArr(const int j_rod, const TRod &rod) {
    const double rLen = rod.length; // Vector of rod
    const double *rUVec = rod.direction;
    const double *rCenter = rod.pos;
    double rVec[3], // Rod length vector
        rMinus[3],  // Rod minus end position vector
        rPlus[3],   // Rod plus end position vector
        sepVec[3];  // Vector of rod center to protein
    for (int i = 0; i < 3; ++i) {
        rVec[i] = rLen * rUVec[i];
        rMinus[i] = rCenter[i] - (.5 * rVec[i]);
        rPlus[i] = rCenter[i] + (.5 * rVec[i]);
        sepVec[i] = pos_[i] - rCenter[i];
    }
    // the perpendicular distance & position from protein ends to rod
    double pointMin[3];
    distMinArr_[j_rod] = dist_point_seg(pos_, rMinus, rPlus, pointMin);
    // Closest point of end_pos along rod axis from rod center.
    double mu0 = dot3(sepVec, rUVec);
    muArr_[j_rod] = mu0;
    // Perpendicular distance away from rod axis
    distPerpArr_[j_rod] = sqrt(dot3(sepVec, sepVec) - SQR(mu0));
}

/*! \brief Calculate the probability of a head to bind to surrounding rods.
 *
 * \param rods Array of rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param bindFactor Binding factor of head to rods \return
 * void, Changes prob_tot_ variable of KMC this object
 */
template <typename TRod>
void KMC<TRod>::CalcTotProbsUS(const TRod *const *rods,
                               const std::vector<int> &uniqueFlagJ,
                               const std::vector<double> &bindFactors) {
    // Make sure that KMC was properly initialized.
    // bindFactorsArr_ = bindFactors;
    assert(r_cutoff_ > 0);
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < rods_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0) {
            rods_probs_[j_rod] =
                CalcProbUS(j_rod, *(rods[j_rod]), bindFactors[j_rod]);
            prob_tot_ += rods_probs_[j_rod];
        }
    }
}

/*! \brief Calculate the probability of a head to bind to one rod.
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::CalcProbUS(const int j_rod, const TRod &rod,
                             const double bindFactor) {
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
    double prefact = bindFactor * dt_ / (4.0 / 3.0 * M_PI * CUBE(r_cutoff_));

    // TODO Make sure max binding factor is less than 1
    // check(prefact*2*rcutoff < 1.)

    UpdateRodDistArr(j_rod, rod);
    double distMinSQR = SQR(distMinArr_[j_rod]);
    double r_cutSQR = SQR(r_cutoff_);

    // Find length of line that goes through binding radius
    if (distMinSQR < r_cutSQR) {
        // Closest point of end_pos along rod from rod center.
        double distPerpSQR = SQR(distPerpArr_[j_rod]);
        double tLen = rod.length;
        double mueff = ABS(muArr_[j_rod]);
        double a = sqrt(r_cutSQR - distPerpSQR);
        double max = mueff + a;
        if (max > 0.5 * tLen)
            max = 0.5 * tLen;
        else if (max < -0.5 * tLen)
            max = -0.5 * tLen;

        double min = mueff - a;
        if (min > 0.5 * tLen)
            min = 0.5 * tLen;
        else if (min < -0.5 * tLen)
            min = -0.5 * tLen;

        return prefact * (max - min);
    } else {
        return 0;
    }
}

/*! \brief Calculate the probability of a head unbinding from a rod it
 * is attached too.
 *
 *  Head must be attached to a rod.
 *
 * \param unbindFactor
 * \return void, Changes the prob_tot_ variable of this object.
 */
template <typename TRod>
void KMC<TRod>::CalcProbSU(const double unbindFactor) {
    prob_tot_ = unbindFactor * dt_;
}

/*! \brief Calculate the total probability of an unbound head to bind to
 * surrounding rods when other head is attached to another rod.
 *
 *  One head must be bound and its position must be stored in this pos_
 * variable.
 *
 * \param rods Array of surrounding rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param k_spring Spring constant between connected ends
 * \param eqLen Equilibrium length of spring connecting ends
 * \param bindFactor Binding factor of head to rods
 * \return void, Changes prob_tot_ variable of this object
 */
template <typename TRod>
void KMC<TRod>::CalcTotProbsSD(const TRod *const *rods,
                               const std::vector<int> &uniqueFlagJ,
                               const int boundID, const double lambda,
                               const double kappa, const double beta,
                               const double restLen,
                               const std::vector<double> &bindFactors) {
    // Make sure that KMC was properly initialized.
    assert(r_cutoff_ > 0);
    // Stash bind factors for position calculation later
    // bindFactorsArr_ = bindFactors;
    // Sum total probability of binding to surrounding rods
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < rods_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0 && rods[j_rod]->gid != boundID) {
            if (LUTablePtr_) {
                rods_probs_[j_rod] =
                    LUCalcProbSD(j_rod, *(rods[j_rod]), bindFactors[j_rod]);
            } else {
                rods_probs_[j_rod] =
                    CalcProbSD(j_rod, *(rods[j_rod]), lambda, kappa, beta,
                               restLen, bindFactors[j_rod]);
            }
            prob_tot_ += rods_probs_[j_rod];
        }
    }
    return;
}

/*! \brief Calculate the probability of a head to bind to one rod when one
 * head is already bound.
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param lambda Load sensitivity of unbinding
 * \param kappa Spring constant of protein when stretched
 * \param beta Inverse product of Boltzmann's constant and temperature
 * \param restLen Length of protein when neither compressed or stretched
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::CalcProbSD(const int j_rod, const TRod &rod,
                             const double lambda, const double kappa,
                             const double beta, const double restLen,
                             const double bindFactor) {
    // // Find and add shortest distance to DistPerp array and the associated
    // // locations along rod.
    UpdateRodDistArr(j_rod, rod);

    double distMinSQR = SQR(distMinArr_[j_rod]);
    double r_cutSQR = SQR(r_cutoff_);
    const double D = 2. * rod.radius;

    double result;
    if (r_cutSQR < distMinSQR)
        result = 0;
    else {
        // Non-dimensionalize all exponential factors
        const double distPerp =
            distPerpArr_[j_rod] / D; // Perpendicular distance to rod axis
        const double mu0 = muArr_[j_rod] / D;
        const double M = (1 - lambda) * .5 * kappa * beta * D * D;
        const double ell = restLen / D;
        // Limits of integration
        const double lim0 = -mu0 - 0.5 * (rod.length / D);
        const double lim1 = -mu0 + 0.5 * (rod.length / D);
        lims_[j_rod].first = lim0;
        lims_[j_rod].second = lim1;
        // It's integrating time!
        result = integral(distPerp, lim0, lim1, M, ell);
    }
    return bindFactor * result * dt_;
}

/*! \brief Calculate the probability of a head to bind to one rod when
 * one head is already bound.using a pre-made lookup table
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param lambda Load sensitivity of unbinding
 * \param kappa Spring constant of protein when stretched
 * \param beta Inverse product of Boltzmann's constant and temperature
 * \param restLen Length of protein when neither compressed or stretched
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::LUCalcProbSD(const int j_rod, const TRod &rod,
                               const double bindFactor) {
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table not initialized ***"
                  << std::endl;
        exit(1);
    }
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
    UpdateRodDistArr(j_rod, rod);

    double result;
    // Bypass probability calculation if protein is too far away
    if (SQR(r_cutoff_) < SQR(distMinArr_[j_rod]))
        result = 0;
    else {
        double mu0 = muArr_[j_rod];
        double distPerp = distPerpArr_[j_rod];   // Perpendicular distance to
                                                 // rod, changes if croslinker
                                                 // is past the rod end point.
        double lim0 = -mu0 - (0.5 * rod.length); // Lower limit of integral
        double lim1 = -mu0 + (0.5 * rod.length); // Upper limit of integral
        lims_[j_rod].first = lim0;
        lims_[j_rod].second = lim1;
        // Get value of integral from end to first bound site.
        //  Since lookup table function is odd you can simply negate the
        //  negative limits.
        // WARNING: LUT stores only half of the table [0,sMax] by symmetry
        // Only the difference between two bounds are needed here.
        // both shifted by a constant
        double term0 = LUTablePtr_->Lookup(distPerp, fabs(lim0)) *
                       ((lim0 < 0) ? -1.0 : 1.0);
        double term1 = LUTablePtr_->Lookup(distPerp, fabs(lim1)) *
                       ((lim1 < 0) ? -1.0 : 1.0);
        result = (term1 - term0);
    }
    return bindFactor * result * dt_;
}

/*! \brief Calculate probability of head unbinding from rod
 *
 *
 * \param unbindFactor Unbinding factor for head that is unbinding
 * \return void, Changes the prob_tot_ membr
 */
template <typename TRod>
void KMC<TRod>::CalcProbDS(const double unbindFactor) {
    prob_tot_ = unbindFactor * dt_;
    return;
}

/*! \brief Find which MT head will bind to from stored probabilities
 *
 * \param rods Array of rod object pointers
 * \param &bindPos Location along rod relative to rod center that head
 * binds
 * \param roll Uniformally generated random number from 0 to 1
 *
 * \return ID of bound MT, changes bindPos to reflect location of bound MT
 */
template <typename TRod>
int KMC<TRod>::whichRodBindUS(const TRod *const *rods, double &bindPos,
                              double roll) {
    assert(roll <= 1.0 && roll >= 0.);
    // Rescale to total binding probability range
    roll *= prob_tot_;

    // assert probabilities are not zero
    double pos_roll = 0.0;

    int i = 0;
    for (auto prob : rods_probs_) {
        if ((pos_roll + prob) > roll) {
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else {
            pos_roll += prob;
            i++;
        }
    }
    if (i == rods_probs_.size()) {
        // Roll given was too large, CHECK KMC step functions
        return -1;
    }
    const double half_length = .5 * rods[i]->length; // half length of rod
    const double mu = muArr_[i]; // Closest position of particle along rod
                                 // from rod center
    const double distMin = distMinArr_[i];   // Distance from rod
    const double distPerp = distPerpArr_[i]; // Distance from rod
    double bind_range = sqrt(SQR(r_cutoff_) - SQR(distPerp));
    if (std::isnan(bind_range))
        bind_range = 0.0;

    // double diff = sqrt(SQR(distMin) - SQR(distPerp));
    double range_m,
        range_p; // Minus and plus binding limits with respect to rod center
    if ((ABS(mu) - (half_length)) > SMALL) { // Protein is beyond rod
        if (mu < 0) {                        // Protein is beyond the minus end
            range_m = -half_length;
            range_p = ((mu + bind_range) > half_length) ? half_length
                                                        : (mu + bind_range);
        } else { // Protein is beyond the plus end
            range_m = ((mu - bind_range) > half_length) ? half_length
                                                        : (mu - bind_range);
            range_p = half_length;
        }
    } else {
        range_m = ((mu - bind_range) < -half_length) ? -half_length
                                                     : (mu - bind_range);
        // Closest particle can bind to rod plus end
        range_p =
            ((mu + bind_range) > half_length) ? half_length : (mu + bind_range);
    }
    // Find location on rod between plus and minus range
    bindPos = (pos_roll * (range_p - range_m)) + range_m;
    assert(bindPos >= -(half_length + 10e-8) &&
           bindPos <= (half_length + 10e-8));
    return i;
}

/*! \brief Converts 3 uniformaly chosen random numbers into spherical
 *   coordinates where head will be released relative to original bound
 * position.
 *
 * \param R Radius of sphere that protein will unbind to
 * \param diffConst Diffusion constant
 * \param rollVec Uniformally generated random number 3 vector all from 0
 * to 1.
 * \return Return new position vector of head once detatched.
 */
template <typename TRod>
void KMC<TRod>::whereUnbindSU(double R, double diffConst, double rollVec[3],
                              double pos[3]) {
    assert(rollVec[0] >= 0 && rollVec[0] <= 1.0);
    assert(rollVec[1] >= 0 && rollVec[1] <= 1.0);
    assert(rollVec[2] >= 0 && rollVec[2] <= 1.0);

    double avg_dist = getDiffRadius(diffConst);
    r_cutoff_ = avg_dist > R ? avg_dist : R;

    double r = R * std::cbrt(rollVec[0]);
    double costheta = 2. * rollVec[1] - 1.;
    double theta = acos(costheta);
    double phi = 2. * M_PI * rollVec[2];
    pos[0] = r * sin(theta) * cos(phi) + pos_[0];
    pos[1] = r * sin(theta) * sin(phi) + pos_[1];
    pos[2] = r * costheta + pos_[2];
    return;
}

/*! \brief Bind unbound protein head to rod close by.
 *
 *
 * \param roll Uniform random number between 0 and 1.
 * \return Index number of rod that head was bound too.
 */
template <typename TRod>
int KMC<TRod>::whichRodBindSD(double &bindPos, double roll) {
    assert(roll >= 0. && roll <= 1.0);
    roll *= prob_tot_; // Rescale to non-normalized value
    double pos_roll = 0.0;
    int i = 0; // Index of rods
    for (auto prob : rods_probs_) {
        if ((pos_roll + prob) > roll) { // Found rod to bind to
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else { // Keep searching for binding rod
            pos_roll += prob;
            i++;
        }
    }
    if (i == rods_probs_.size()) {
        return -1;
    }
    bindPos = RandomBindPosSD(i, pos_roll);
    return i;
}

/*! \brief Where on a rod will a protein bind if already rod to another rod
 *
 * \param j_rod Parameter description
 * \param roll A uniformly generated random number between 0-1
 * \return Return parameter description
 */
template <typename TRod>
double KMC<TRod>::RandomBindPosSD(int j_rod, double roll) {
    assert(roll <= 1.0 && roll >= 0.);
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table pointer not initialized ***"
                  << std::endl;
        exit(1);
    }
    // Limits set at the end of the rod relative to closest point of bound head
    // to unbound rod (mu)
    double lim0 = lims_[j_rod].first;
    double lim1 = lims_[j_rod].second;

    // Lookup table parameter: perpendicular distance away from rod
    const double distPerp = distPerpArr_[j_rod];

    // Range of CDF based on position limits where protein can bind
    double pLim0 =
        LUTablePtr_->Lookup(distPerp, fabs(lim0)) * ((lim0 < 0) ? -1.0 : 1.0);
    double pLim1 =
        LUTablePtr_->Lookup(distPerp, fabs(lim1)) * ((lim1 < 0) ? -1.0 : 1.0);

    // Scale random number to be within probability limits
    //      Since lookup table only stores positive values, take the absolute
    //      value of the position limit and then negate probability if the
    //      position limit was less than 0.
    roll = roll * (pLim1 - pLim0) + pLim0;
    double preFact = roll < 0 ? -1. : 1.;
    roll *= preFact; // entry in lookup table must be positive

    double bind_pos = preFact * (LUTablePtr_->ReverseLookup(distPerp, roll));
#ifndef NDEBUG
    printf("roll = %f, expected range = (%f, %f) \n", preFact * roll, pLim0,
           pLim1);
    printf("bind_pos = %f, expected range = (%f, %f) \n", bind_pos, lim0, lim1);
#endif

    // FIXME Reverse lookup table is not perfect. Bind to end of MT if limit is
    // exceeded.
    if (lim0 > bind_pos) {
#ifndef NDEBUG
        printf("Warning: binding position %g exceeds rod limit of %g. Binding "
               "to end \n",
               bind_pos, lim0);
#endif
        bind_pos = lim0;
    } else if (lim1 < bind_pos) {
#ifndef NDEBUG
        printf("Warning: binding position %g exceeds rod limit of %g. Binding "
               "to end \n",
               bind_pos, lim1);
#endif
        bind_pos = lim1;
    }
    // assert(lim0 < bind_pos && lim1 > bind_pos);
    // Translate bound position to be relative to center of rod
    return bind_pos + muArr_[j_rod];
}

/*! \brief Where should head be located after unbinding while other head is
 *  is still bound.
 *
 *  At the moment, head just goes to back to location of other bound head.
 * This should not change the behavior of the KMC step presently. Might be
 * changed in the future when crosslinkers are longer.
 *
 * \param boundPos[] Parameter description
 * \param pos[] Parameter description
 * \return Return parameter description
 */
template <typename TRod>
void KMC<TRod>::whereUnbindDS(const double boundPos[], double pos[]) {
    for (int i = 0; i < 3; ++i) {
        pos[i] = boundPos[i];
    }
    return;
};

#endif /* KMC_HPP */
