/**
 * @author      : alamson (alamson@ccblin059.flatironinstitute.org)
 * @file        : kmc_overload
 * @created     : Saturday Mar 06, 2021 11:31:21 EST
 */

#ifndef KMC_OVERLOAD_HPP

#define KMC_OVERLOAD_HPP

// Overload with NpjSphere=0 (backwards compatibility)
KMC(const double *pos, const int NpjRod, const double r_cutoff,
    const double diffConst, const double dt)
    : dt_(dt), LUTablePtr_(nullptr) {
    setPos(pos);

    // Find average diffusion distance and compare to given r_cutoff
    double avg_dist = getDiffRadius(diffConst);
    r_cutoff_ = (avg_dist > r_cutoff) ? avg_dist : r_cutoff;
    bind_vol_ = M_PI * CUBE(r_cutoff_) / .75;

    // Set size of rod-sized vectors
    rod_probs_.resize(NpjRod, 0);
    distMinArrRod_.resize(NpjRod, 0);
    distPerpArr_.resize(NpjRod, 0);
    muArr_.resize(NpjRod, 0);
    lims_.resize(NpjRod);

    assert(r_cutoff_ > 0);
}

// Overload with NpjSphere=0 (backwards compatibility)
KMC(const double *pos, const int NpjRod, const double dt,
    const LookupTable *LUTablePtr)
    : dt_(dt), LUTablePtr_(LUTablePtr) {
    setPos(pos);
    r_cutoff_ = LUTablePtr_->getLUCutoff();
    bind_vol_ = LUTablePtr_->getBindVolume();

    // Set size of rod-sized vectors
    rod_probs_.resize(NpjRod, 0);
    distMinArrRod_.resize(NpjRod, 0);
    distPerpArr_.resize(NpjRod, 0);
    muArr_.resize(NpjRod, 0);
    lims_.resize(NpjRod);
}

// Overload with NpjSphere=0 (backwards compatibility)
KMC(const double *pos, const int NpjRod, const double r_cutoff, const double dt)
    : r_cutoff_(r_cutoff), dt_(dt), LUTablePtr_(nullptr) {
    setPos(pos);

    // Set size of rod-sized vectors
    rod_probs_.resize(NpjRod, 0);
    distMinArrRod_.resize(NpjRod, 0);
    distPerpArr_.resize(NpjRod, 0);
    muArr_.resize(NpjRod, 0);
    lims_.resize(NpjRod);
}

// Overload for empty sphere vector
void CalcTotProbsUS(const std::vector<const TRod *> &rods,
                    const std::vector<double> &bindFactors) {
    std::vector<const TSphere *> spheres;
    CalcTotProbsUS(rods, spheres, bindFactors);
}

// Overload with empty spheres vector
void LUCalcTotProbsSD(const std::vector<const TRod *> &rods, const int boundID,
                      const std::vector<double> &bindFactors) {
    std::vector<const TSphere *> spheres;
    LUCalcTotProbsSD(rods, spheres, boundID, bindFactors);
}

// Overload with array (backwards compatibility)
// int whichRodBindUS(TRod** rods, double &bindPos, double roll) {
int whichRodBindUS(const TRod *const *rods, double &bindPos, double roll) {
    // Convert rods to vector
    std::vector<const TRod *> rodsVec(rods, rods + rod_probs_.size());
    return whichRodBindUS(rodsVec, bindPos, roll);
}

#endif /* end of include guard KMC_OVERLOAD_HPP */

