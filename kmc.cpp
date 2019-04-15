/**
 * @file kmc.cpp
 * @brief KMC object and speciallized call functions to implement
 * two stage kmc implementation with code base.
 */
#include "kmc.hpp"

#include <cmath>

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
void KMC::UpdateRodDistArr(const int j_rod, const SylinderNearEP &rod) {
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
    // PS::F64vec3 muVec = pointMin - rCenter;
    // Closest point of end_pos along rod axis from rod center.
    double mu0 = dot3(sepVec, rUVec);
    muArr_[j_rod] = mu0;
    // Perpendicular distance away from rod axis
    distPerpArr_[j_rod] = sqrt(dot3(sepVec, sepVec) - SQR(mu0));

    assert((distPerpArr_[j_rod] - distMinArr_[j_rod]) > 10e-8);
}

/*! \brief Calculate the probability of a head to bind to surrounding rods.
 *
 * \param ep_j Array of rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param bindFactor Binding factor of head to rods \return
 * void, Changes prob_tot_ variable of KMC this object
 */
void KMC::CalcTotProbsUS(const SylinderNearEP *const *ep_j,
                         const std::vector<int> &uniqueFlagJ,
                         const double bindFactor) {
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < ep_j_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0) {
            ep_j_probs_[j_rod] = CalcProbUS(j_rod, *(ep_j[j_rod]), bindFactor);
            prob_tot_ += ep_j_probs_[j_rod];
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
double KMC::CalcProbUS(const int j_rod, const SylinderNearEP &rod,
                       const double bindFactor) {
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
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

        return bindFactor * (max - min) / (4.0 / 3.0 * M_PI * CUBE(r_cutoff_));
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
void KMC::CalcProbSU(const double unbindFactor) { prob_tot_ = unbindFactor; }

/*! \brief Calculate the total probability of an unbound head to bind to
 * surrounding rods when other head is attached to another rod.
 *
 *  One head must be bound and its position must be stored in this pos_
 * variable.
 *
 * \param ep_j Array of surrounding rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param k_spring Spring constant between connected heads
 * \param eqLen Equilibrium length of spring connecting heads
 * \param bindFactor Binding factor of head to rods
 * \return void, Changes tot_prob_ variable of this object
 */
void KMC::CalcTotProbsSD(const SylinderNearEP *const *ep_j,
                         const std::vector<int> &uniqueFlagJ, const int boundID,
                         const double lambda, const double kappa,
                         const double beta, const double restLen,
                         const double bindFactor) {
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < ep_j_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0 && ep_j[j_rod]->gid != boundID) {
            if (LUTablePtr_) {
                // ep_j_probs_[j_rod] = LUCalcProbSD(j_rod, *(ep_j[j_rod]),
                // kappa, eqLen, bindFactor);
                ep_j_probs_[j_rod] =
                    LUCalcProbSD(j_rod, *(ep_j[j_rod]), bindFactor);
            } else {
                ep_j_probs_[j_rod] =
                    CalcProbSD(j_rod, *(ep_j[j_rod]), lambda, kappa, beta,
                               restLen, bindFactor);
            }
            prob_tot_ += ep_j_probs_[j_rod];
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
double KMC::CalcProbSD(const int j_rod, const SylinderNearEP &rod,
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
    return bindFactor * result;
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
double KMC::LUCalcProbSD(const int j_rod, const SylinderNearEP &rod,
                         const double bindFactor) {
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table not initialized ***"
                  << std::endl;
        exit(1);
    }
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
    UpdateRodDistArr(j_rod, rod);

    double D = 2. * rod.radius;
    double result;
    // Bypass probability calculation if protein is too far away
    if (SQR(r_cutoff_) < SQR(distMinArr_[j_rod]))
        result = 0;
    else {
        double mu0 = muArr_[j_rod] / D;
        double distPerp = distPerpArr_[j_rod] / D; // Perpendicular distance to
                                                   // rod, changes if croslinker
                                                   // is past the rod end point.
        double lim0 = -mu0 - (0.5 * rod.length / D); // Lower limit of integral
        double lim1 = -mu0 + (0.5 * rod.length / D); // Upper limit of integral
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
        // Integral must have dimensions of length
        result = (term1 - term0) * D;
    }
    return bindFactor * result;
}

/*! \brief Calculate probability of head unbinding from rod
 *
 *
 * \param unbindFactor Unbinding factor for head that is unbinding
 * \return void, Changes the prob_tot_ membr
 */
void KMC::CalcProbDS(const double unbindFactor) {
    prob_tot_ = unbindFactor;
    return;
}

/*! \brief Find which MT head will bind to from stored probabilities
 *
 * \param ep_j Array of rod object pointers
 * \param &bindPos Location along rod relative to rod center that head
 * binds
 * \param roll Uniformally generated random number from 0 to 1 \return
 * Return ID of bound MT
 */
int KMC::whichRodBindUS(const SylinderNearEP *const *ep_j, double &bindPos,
                        double roll) {
    // assert probabilities are not zero
    double pos_roll = 0.0;

    int i = 0;
    for (auto prob : ep_j_probs_) {
        if ((pos_roll + prob) > roll) {
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else {
            pos_roll += prob;
            i++;
        }
    }
    if (i == ep_j_probs_.size()) {
        // Roll given was too large, CHECK KMC step functions
        return -1;
    }
    const double half_length = .5 * ep_j[i]->length; // half length of rod
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
 * \param rollVec Uniformally generated random number 3 vector all from 0
 * to 1. \return Return new position vector of head once detatched.
 */
void KMC::whereUnbindSU(double R, double rollVec[3], double pos[3]) {
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
int KMC::whichRodBindSD(double &bindPos, double roll) {
    double pos_roll = 0.0;
    int i = 0; // Index of rods
    for (auto prob : ep_j_probs_) {
        if ((pos_roll + prob) > roll) { // Found rod to bind to
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else { // Keep searching for binding rod
            pos_roll += prob;
            i++;
        }
    }
    if (i == ep_j_probs_.size()) {
        return -1;
    }
    bindPos = RandomBindPosSD(i, pos_roll);
    return i;
}

/*! \brief Where on a rod will a protein bind if already rod to another rod
 *
 * \param j_rod Parameter description
 * \param roll A uniformaly generated random number between 0-1
 * \return Return parameter description
 */
double KMC::RandomBindPosSD(int j_rod, double roll) {
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table pointer not initialized ***"
                  << std::endl;
        exit(1);
    }
    const double D = LUTablePtr_->getTubuleDiameter();
    double lim0 = lims_[j_rod].first;
    double lim1 = lims_[j_rod].second;
    // TODO: Clear up this
    // Lookup table parameter
    const double distPerp = distPerpArr_[j_rod] / D;
    // Range of CDF based on position limits where protein must bind
    double pLim0 =
        LUTablePtr_->Lookup(distPerp, fabs(lim0)) * ((lim0 < 0) ? -1.0 : 1.0);
    double pLim1 =
        LUTablePtr_->Lookup(distPerp, fabs(lim1)) * ((lim1 < 0) ? -1.0 : 1.0);
    // Scale random number to be within pLims
    roll = roll * (pLim1 - pLim0) + pLim0;
    double sbound = 0; // Reset x position parameter
    double preFact = roll < 0 ? -1. : 1.;
    roll *= preFact; // entry in lookup table must be positive
    return preFact * (D * LUTablePtr_->ReverseLookup(distPerp, roll)) +
           muArr_[j_rod];
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
void KMC::whereUnbindDS(const double boundPos[], double pos[]) {
    for (int i = 0; i < 3; ++i) {
        pos[i] = boundPos[i];
    }
    return;
}

/**
 * @brief: Perform kinetic monte carlo step of protein object when no heads
 * are attached.
 *
 * @param: &pData Refernce for head object
 * @param: Npj Number of rods near head object
 * @param: ep_j Array of rod object pointers near head
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 *
 * @return: void, Change protein data if head binds
 */
void KMC_U(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double roll,
           ProteinBindStatus &pBind) {
    // Assert the heads are not attached
    assert(pData.getBindID(0) == ID_UB && pData.getBindID(1) == ID_UB);
    // Create KMC objects for probability and step calculations
    double rcut01 = pData.getRcutUS();
    // Create KMC objects
    KMC kmc_end0(pData.bind.pos, Npj, rcut01);
    KMC kmc_end1(pData.bind.pos, Npj, rcut01);
    // Loop over rods and calculate binding probabilities
    kmc_end0.CalcTotProbsUS(ep_j, uniqueFlagJ, pData.getBindingFactorUS(0, dt));
    kmc_end1.CalcTotProbsUS(ep_j, uniqueFlagJ, pData.getBindingFactorUS(1, dt));
    double totBindProb = kmc_end0.getTotProb() + kmc_end1.getTotProb();

    int head_bound = -1; // Head that gets bound, -1 = No binding occurs
    if (totBindProb > 1.0) {
        head_bound = (roll < (kmc_end0.getTotProb() / totBindProb)) ? 0 : 1;
        std::cerr << " *** Warning: Probability of head binding is >1 ("
                  << totBindProb << "). Change time step to prevent this! ***"
                  << std::endl;
    } else if (roll < totBindProb) {
        head_bound = (roll < kmc_end0.getTotProb()) ? 0 : 1;
    } else { // No binding occured
#ifndef NDEBUG
        if (totBindProb > 0)
            printf("Total binding probability: %f \n", totBindProb);
#endif
        return;
    }
    assert(head_bound != -1);
    int rj; // Sylinder j;
    double bindPos = 0;
    if (head_bound == 0) {
        rj = kmc_end0.whichRodBindUS(ep_j, bindPos, roll);
    } else {
        // Shift and scale random roll so it neglects end 0 probabilities
        roll = roll - kmc_end0.getTotProb();
        rj = kmc_end1.whichRodBindUS(ep_j, bindPos, roll);
    }
    // Should not give -1 since roll is within the total binding probability
    // If failure, make sure rolls are shifted properly
    assert(rj != -1);
    // Bind head to rod
    auto &rod = *(ep_j[rj]);
    PS::F64vec3 rVec = rod.getPos();
    double rPos[3] = {rVec[0], rVec[1], rVec[2]};
#ifndef NDEBUG
    printf("U->S Binding\n");
#endif

    pBind.setBind(head_bound, rod.gid, rod.direction, rPos, bindPos, rod.length,
                  rod.rank);
    return;
}

/**
 * @brief: Perform kinetic monte carlo step of protein with 1 head attached.
 *
 * @param: &pData Refernce for protein object
 * @param: Npj Number of rods near protein object
 * @param: ep_j Array of rod objects near protein object
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 * @param: rollVec
 *
 * @return: void, Change protein data if it binds or unbinds
 */
void KMC_S(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
           double rollVec[3], ProteinBindStatus &pBind) {
    double roll = rollVec[0];
    // Find out which head is bound
    int head_bound = (pData.getBindID(0) == ID_UB) ? 1 : 0;

    const ProteinType *pType = &pData.property;

    // Set up KMC objects and calculate probabilities
    KMC kmc_unbind(pData.bind.posEndBind[head_bound]);
    KMC kmc_bind(pData.bind.posEndBind[head_bound], Npj, pData.getRcutSD(),
                 pType->LUTablePtr);
    kmc_unbind.CalcProbSU(pData.getUnbindingFactorSU(head_bound, dt));
    kmc_bind.CalcTotProbsSD(ep_j, uniqueFlagJ, pData.getBindID(head_bound),
                            pType->lambda, pType->kappa, 1. / KBT,
                            pType->freeLength,
                            pData.getBindingFactorSD(head_bound, dt));
    // Get total probability of changing protein state
    double totProb = kmc_bind.getTotProb() + kmc_unbind.getTotProb();

    int head_activate = -1; // No head activated
    if (totProb > 1.0) {    // Probability of KMC is greater than one, normalize
        head_activate = (roll < (kmc_unbind.getTotProb() / totProb))
                            ? head_bound
                            : 1 - head_bound;
        std::cerr << " *** Warning: Probability of head binding or unbinding "
                  << "is >1 (" << totProb << ")."
                  << "Change time step to prevent this! ***" << std::endl;
    } else if (roll <
               totProb) { // Choose which action to perform, bind or unbind
        head_activate =
            (roll < kmc_unbind.getTotProb()) ? head_bound : 1 - head_bound;
    } else { // No head binds or unbinds
        return;
    }
    assert(head_activate != -1);
    // Change status of activated head
    if (head_bound == head_activate) { // Bind unbound head
        rollVec[0] = roll / kmc_unbind.getTotProb();
        double pos[3] = {};
        kmc_unbind.whereUnbindSU(pData.getRcutUS(), rollVec, pos);
        pBind.setUnBindPos(pos);
    } else if (head_activate == 1 - head_bound) { // Unbind bound head
        roll = roll - kmc_unbind.getTotProb();
        double bindPos; // Position on rod where protein will bind,
                        // passed by reference.
        // Pick rod to bind to
        // Should not give -1 since roll is within the total binding
        // probability If failure, make sure rolls are shifted properly
        int rj = kmc_bind.whichRodBindSD(bindPos, roll);
        assert(rj != -1);
        const auto &rod = *(ep_j[rj]);
        // Get position of rod
        PS::F64vec3 rVec = rod.getPos();
        double rPos[3] = {rVec[0], rVec[1], rVec[2]};
#ifndef NDEBUG
        printf("S->D Binding\n");
#endif
        // Bind protein to rod
        pBind.setBind(head_activate, rod.gid, rod.direction, rPos, bindPos,
                      rod.length, rod.rank);
    }
    return;
}

/**
 * @brief: Perform kinetic monte carlo step of protein with 2 heads of
 * protein object attached.
 *
 * @param: &pData Refernce for protein object
 * @param: Npj Number of rods near protein object
 * @param: ep_j Array of rod objects near protein object
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 *
 * @return: void, Change protein data if protein unbinds
 */
void KMC_D(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
           double roll, ProteinBindStatus &pBind) {
    KMC kmc_UB0(pData.bind.posEndBind[0]);
    KMC kmc_UB1(pData.bind.posEndBind[1]);
    kmc_UB0.CalcProbDS(pData.getUnbindingFactorDS(0, dt, KBT));
    kmc_UB1.CalcProbDS(pData.getUnbindingFactorDS(1, dt, KBT));
    double totProb = kmc_UB0.getTotProb() + kmc_UB1.getTotProb();

    int head_activate = -1; // No head activated
    if (totProb > 1.0) {    // Probability of KMC is greater than one, normalize
        std::cerr << " *** Warning: Probability of head binding or unbinding "
                  << "is >1 (" << totProb << ")."
                  << "Change time step to prevent this! ***" << std::endl;
        head_activate = (roll < (kmc_UB0.getTotProb() / totProb)) ? 0 : 1;

    } else if (roll < totProb) { // Choose which head to unbind
        head_activate = (roll < kmc_UB0.getTotProb()) ? 0 : 1;

    } else { // No head unbinds
        return;
    }
    assert(head_activate != -1);
    pBind.setUnBind(head_activate);

    return;
}
