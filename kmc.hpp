#ifndef KMC_HPP
#define KMC_HPP

#include "helpers.hpp"
#include "lookup_table.hpp"
#include "macros.hpp"

#include "Protein/ProteinBindStatus.hpp"
#include "Protein/ProteinData.hpp"
#include "Protein/ProteinType.hpp"
#include "SylinderNear.hpp"

#include <array>
#include <cassert>

// LookupTable BuildLookupTable(const double kMT, const double lenMP0,
// const double y_cutoff);

//// KMC step for unbound protein
// void KMC_U(const ProteinData &pData, const int Npj,
//           const TRod *const *rods,
//           const std::vector<int> &uniqueFlagJ, double dt, double roll,
//           ProteinBindStatus &pBind);
//// KMC step for singly bound to either unbound or doubly bound.
// void KMC_S(const ProteinData &pData, const int Npj,
//           const TRod *const *rods,
//           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
//           double rollVec[], ProteinBindStatus &pBind);
//// KMC step for doubly
// void KMC_D(const ProteinData &bindData, const int Npj,
//           const TRod *const *rods,
//           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
//           double roll, ProteinBindStatus &pBind);

template <class TRod>
class KMC {
  private:
    // Probabilities
    double prob_tot_;
    double r_cutoff_;
    std::vector<double> rods_probs_; ///< binding probabilities

    // Spatial variables
    double pos_[3];                   ///< position of motor head
    std::vector<double> distMinArr_;  ///< min dist to rod segment
    std::vector<double> distPerpArr_; ///< min (perp) dist to rod line
    std::vector<double> muArr_;
    ///< Closest points of proteins to rods along the rods axis with
    ///< respect to rod center and direction.

    std::vector<std::pair<double, double>> lims_;
    ///< Bounds of binding positions on rods with respect to closest
    ///< point of rod axis. Used in lookup table value acquisition.

    LookupTable *LUTablePtr_ = nullptr; ///< lookup table

  public:
    /******************
     *  Constructors  *
     ******************/

    KMC() {}

    KMC(const double *pos) { setPos(pos); }

    KMC(const double *pos, const int Npj, const double r_cutoff) {
        setPos(pos);
        r_cutoff_ = r_cutoff;
        prob_tot_ = 0;
        rods_probs_.resize(Npj, 0);
        distMinArr_.resize(Npj, 0);
        distPerpArr_.resize(Npj, 0);
        muArr_.resize(Npj, 0);
        lims_.resize(Npj);
    }

    KMC(const double *pos, const int Npj, const double r_cutoff,
        LookupTable *LUTablePtr) {
        setPos(pos);
        r_cutoff_ = r_cutoff;
        prob_tot_ = 0;
        rods_probs_.resize(Npj, 0);
        distMinArr_.resize(Npj, 0);
        distPerpArr_.resize(Npj, 0);
        muArr_.resize(Npj, 0);
        lims_.resize(Npj);
        LUTablePtr_ = LUTablePtr;
    }

    /*************************************
     *  Calculate probability functions  *
     *************************************/

    void CalcTotProbsUS(const TRod *const *rods,
                        const std::vector<int> &uniqueFlagJ,
                        const double bindFactor);

    double CalcProbUS(const int j_bond, const TRod &rod,
                      const double bindFactor);

    void CalcProbSU(const double unbindFactor);

    void CalcTotProbsSD(const TRod *const *rods,
                        const std::vector<int> &uniqueFlagJ, const int boundID,
                        const double lambda, const double kappa,
                        const double beta, const double restLen,
                        const double bindFactor);

    double LUCalcProbSD(const int j_rod, const TRod &rod,
                        const double bindFactor);

    double CalcProbSD(const int j_rod, const TRod &rod, const double lambda,
                      const double kappa, const double beta,
                      const double restLen, const double bindFactor);

    void CalcProbDS(const double unbindFactor);

    void UpdateRodDistArr(const int j_bond, const TRod &rod);

    double getTotProb() { return prob_tot_; }

    int whichRodBindUS(const TRod *const *rods, double &bindPos, double roll);

    void whereUnbindSU(double R, double rollVec[3], double pos[3]);

    int whichRodBindSD(double &bindPos, double roll);

    double RandomBindPosSD(int j_bond, double roll);

    void whereUnbindDS(const double boundPos[], double pos[]);

    /*******************
     *  Get functions  *
     *******************/

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

template <class TRod>
TRod MockSylinder(int id);

#endif /* KMC_HPP */
