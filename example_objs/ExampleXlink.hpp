
/**
 * @file ProteinBindStatus.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-02-07
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef EXAMPLEXLINK_HPP_
#define EXAMPLEXLINK_HPP_

#include "KMC/lookup_table.hpp"
#include "KMC/macros.hpp"
#include <cassert>
#include <limits>
#include <type_traits>

constexpr int IDUB = -1; // the default Unbind state
constexpr double NAN_D = std::numeric_limits<double>::quiet_NaN();
constexpr double EPS_D = std::numeric_limits<double>::epsilon();

/**
 * @brief time-varying ProteinBindStatus
 *
 * This class offers only basic setBind() and setUnbind() functions.
 * Other protein behaviors are defined in ProteinData
 *
 */
class ExampleXlink {
  public:
    // TODO: Test this
    double pos[3] = {NAN_D, NAN_D, NAN_D};
    ///< Case1 both unbind -> protein is a point only
    ///< Case2 one bind -> protein is a point only, pos = posEndBind
    ///< Case3 both bind -> pos = mid-point of posEndBind

    // time varying properties
    bool changeBind[2] = {true, true};  ///< flag for if binding status changes
    int idBind[2] = {IDUB, IDUB};       ///< ID of bind MT
    int rankBind[2] = {IDUB, IDUB};     ///< mpi rank of bind MT
    double lenBind[2] = {NAN_D, NAN_D}; ///< length of bind MT
    double distBind[2] = {NAN_D, NAN_D};
    ///< the distance to bind MT center,
    ///< [-lenBind/2,lenBind/2], +towards plus end
    double posEndBind[2][3] = {{NAN_D, NAN_D, NAN_D}, {NAN_D, NAN_D, NAN_D}};
    ///< position (in lab frame) of two ends when bind
    ///< when unbind, the position is NAN_D
    double centerBind[2][3] = {{NAN_D, NAN_D, NAN_D}, {NAN_D, NAN_D, NAN_D}};
    ///< center position of two MTs
    double directionBind[2][3] = {{NAN_D, NAN_D, NAN_D}, {NAN_D, NAN_D, NAN_D}};
    ///< direction of two MTs
    int tag = 0;              ///< user-assigned integer tag for different types
    bool walkOff = true;      ///< walf off the end, i.e., no 'end-dewelling'
    double fixedLocation = 0; ///< in [-1,1]
    double fstall = 1.0;      ///< stall force. unit pN
    double diffUnbound = .1;  ///< Unbounded diffusivity, unit um^2/s
    double lambda = 0;        ///< dimensionless unbinding load sensitivity

    // TODO: Implement these two functions
    bool bindAntiParallel = false; ///< only bind anti parallel MTs
    bool fixedEnd0 = false; ///< end0 is fixed to a MT at the given location

    // per head properties
    double vmax[2]; ///< singly bound max velocity for each end. unit um/s

    // kMC data field
    double freeLength = 0.05; ///< the free length. unit um
    double kappa = 204.7; ///< Spring constant when attached to MT. unit pN/um
    double rc;            ///< the capture radius of protein (user set)
    double eps;           ///< Number of crosslinker binding sites per um of MT
    double ko_s[2];       ///< End turnover rates or KMC_u_s and KMC_s_u steps
    double ko_d[2];       ///< End turnover rates for KMC_s_d and KMC_d_s steps
    double Ka[2];         ///< Association constant when neither head is  bound
    double Ke[2];         ///< Equilibrium constant when other end is bound

    // kMC lookup table
    LookupTable *LUTablePtr = nullptr; ///< Pointer to lookup table for binding

    void echo() const {
        printf("------Protein Type Properties-----\n");
        printf("tag: %d\n", tag);
        printf("walk off: %d\n", walkOff);
        printf("bindAntiParallel: %d\n", bindAntiParallel);
        printf("fixedEnd0: %d\n", fixedEnd0);
        printf("freeLength: %g um\n", freeLength);
        printf("spring lappa: %g pN/um\n", kappa);
        printf("fstall: %g pN\n", fstall);
        printf("capture radius: %g um\n", rc);
        printf("Unbinding load sensitivity: %g \n", lambda);
        printf("Binding site density along MT: %g 1/um \n", eps);
        printf("vmax: %g um/s, %g um/s\n", vmax[0], vmax[1]);
        printf("singly bound turnover rate: %g 1/s,%g 1/s\n", ko_s[0], ko_s[1]);
        printf("doubly bound turnover rate: %g 1/s,%g 1/s\n", ko_d[0], ko_d[1]);
        printf("Ka_s: %g um^3,%g um^3\n", Ka[0], Ka[1]);
        printf("Ke_d dimensionless: %g ,%g \n", Ke[0], Ke[1]);
        printf("----------------------------------\n");
    }

    const double *getPosPtr() const { return pos; }

    /**
     * @brief Set the end to unbind status
     *
     * @param end 0 or 1
     */
    void setUnBind(int end) {
        idBind[end] = IDUB;
        rankBind[end] = IDUB;
        lenBind[end] = NAN_D;
        distBind[end] = NAN_D;
        for (int i = 0; i < 3; i++) {
            directionBind[end][i] = NAN_D;
            centerBind[end][i] = NAN_D;
            posEndBind[end][i] = NAN_D;
        }
    }

    /**
     * @brief Set position for completely unbound protein
     *
     */
    void setUnBindPos(double newPos[3]) {
        pos[0] = newPos[0];
        pos[1] = newPos[1];
        pos[2] = newPos[2];
    }

    /**
     * @brief: setBind status
     *
     * @param: int end
     * @param: gid
     * @param: directionLine
     * @param: centerLine
     * @param: centerDist
     *
     * @return: void
     */
    void setBind(int end, const int gid, const double directionLine[3],
                 const double centerLine[3], const double centerDist,
                 const double length, const int rank) {
        assert(end == 0 || end == 1); // Make sure you choose a viable head
        assert(idBind[end] == IDUB);  // Make sure end is originally unbound
        assert(gid != IDUB);
        idBind[end] = gid;
        distBind[end] = centerDist;
        lenBind[end] = length;
        rankBind[end] = rank;
        std::copy(directionLine, directionLine + 3, directionBind[end]);
        std::copy(centerLine, centerLine + 3, centerBind[end]);
        updatePosEndBind(end);
    }

    /**
     * @brief update posEndBind and pos
     *
     */
    void updateGeometryWithBind() {
        updatePosEndBind(0);
        updatePosEndBind(1);
        updatePosWithEndBind();
    }

    /**
     * @brief calculate the position of end with bind information
     *
     * @param end
     */
    void updatePosEndBind(const int end) {
        if (idBind[end] == IDUB) {
            setUnBind(end);
        } else {
            for (int i = 0; i < 3; i++) { // posEnd = direction * dist + center
                posEndBind[end][i] =
                    directionBind[end][i] * distBind[end] + centerBind[end][i];
            }
        }
    }

    /**
     * @brief
     * Case 1: if doubly bound, set pos to center of two ends.
     * Case 2: if singly bound, set pos to that end.
     * Case 3: if unbound, do nothing
     * This must be called when posEndBind is valid
     */
    void updatePosWithEndBind() {
        if (idBind[0] != IDUB && idBind[1] == IDUB) { // Case 2
            std::copy(posEndBind[0], posEndBind[0] + 3, pos);
        } else if (idBind[0] == IDUB && idBind[1] != IDUB) { // Case 2
            std::copy(posEndBind[1], posEndBind[1] + 3, pos);
        } else if (idBind[0] != IDUB && idBind[1] != IDUB) { // Case 1
            for (int i = 0; i < 3; i++) {
                pos[i] = 0.5 * (posEndBind[0][i] + posEndBind[1][i]);
            }
        }
    }

    /**
     * @brief clamp posEnd to MT range
     *
     * @param end
     */
    void updatePosEndClamp(int end) {
        if (idBind[end] != IDUB) {
            const double lenHalf = lenBind[end] * 0.5;
            if (distBind[end] > lenHalf)
                distBind[end] = lenHalf;
            if (distBind[end] < -lenHalf)
                distBind[end] = -lenHalf;
        }
    }

    double getBindingFactorUS(int e) const {
        assert(e == 0 || e == 1);
        return ko_s[e] * Ka[e] * eps;
    }

    /**
     * @brief Get UnbindingFactor from Singly bound to Unbound
     *
     * @param e
     * @return double
     */
    double getUnbindingFactorSU(int e) const { return ko_s[e]; }

    /**
     * @brief Get BindingFactor from Singly to Doubly bound
     *
     * @param e
     * @return double
     */
    double getBindingFactorSD(int e) const { return ko_d[e] * Ke[e] * eps; }

    /**
     * @brief Get UnbindingFactor from Doubly to Singly bound
     *
     * @param e
     * @return double
     */
    double getUnbindingFactorDS(int e, double KBT) const {
        // if (lambda == 0)
        return ko_d[e];
        // else {
        // double M = lambda * .5 * kappa *
        // SQR(getProteinForceLength() - freeLength) / KBT;
        // return ko_d[e] * dt * exp(M);
        //}
    }

    /**
     * @brief Get cut off calculation range when calculating 0 to 1 head bound
     *
     * @return double
     */
    double getRcutUS() const { return rc; }

    /**
     * @brief Get cut off calculation range when calculating 0 to 1 head bound
     *
     * @return double
     */
    double getDiffU() const { return diffUnbound; }

    /**
     * @brief Get cut off calculation range when calculating 1 to 2 head bound
     * TODO: Test in KMC_Testing
     * @return double
     */
    double getRcutSD() const {
        // consistent tubuleDiameter as the original LUT construction
        double rodDiameter = LUTablePtr->getRodDiameter();
        // use the cutoff length as LUT construction
        return LUTablePtr->getNonDsbound() * rodDiameter;
    }

    /*! \brief Return the ID of the rod protein end is bound to
     *
     *
     * \param end Which protein end
     * \return ID of rod protein is bound to
     */
    int getBindID(int end) const { return idBind[end]; }

    /*! \brief Create a mock xlink for easy testing. Starts unbound at pos
     *  {0, 0, 0}
     *
     * \return void
     */
    void setMockXlink() {
        std::fill(pos, pos + 3, 0);
        Ka[0] = Ka[1] = 1.0;
        Ke[0] = Ke[1] = 1.0;
        ko_s[0] = ko_s[1] = 1.0;
        ko_d[0] = ko_d[1] = 1.0;
        vmax[0] = vmax[1] = 1.0; // max velocity. unit um/s
        kappa = 1.0;             // Spring constant when attached to MT
        diffUnbound = .1;        ///< Unbounded diffusivity, unit um^2/s
        eps = 1.0;
        rc = .5;          // capture radius
        freeLength = 2.0; // the 'bind search' radius
        fstall = 1.0;     ///< stall force. unit pN
        tag = 0;          ///< user-assigned integer tag for different types
        walkOff = true;   ///< walf off the end, i.e., no 'end-dewelling'
        lambda = 0;
    }

    double getExponentFactor() const { return .5 * kappa * (1. - lambda); }
};

static_assert(std::is_trivially_copyable<ExampleXlink>::value, "");
static_assert(std::is_default_constructible<ExampleXlink>::value, "");

#endif
