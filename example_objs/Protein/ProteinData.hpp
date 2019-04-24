/**
 * @file ProteinData.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef PROTEINDATA_HPP_
#define PROTEINDATA_HPP_

#include "ProteinBindStatus.hpp"
#include "ProteinType.hpp"

#include "lookup_table.hpp"
#include "macros.hpp"

//#include "SimToolbox/FDPS/particle_simulator.hpp"
//#include "SimToolbox/Util/EigenDef.hpp"
//#include "SimToolbox/Util/GeoUtil.hpp"
//#include "SimToolbox/Util/IOHelper.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

/**
 * @brief FullParticle type for Protein managed by PS::ParticleSystem
 *
 *  TODO: clear up this class
 * Naming convention:
 * setXXX(), getXXX(): usual set&get
 * updateXXX(): update member state, return void
 * calcXXX() const: calc with some input, return value
 */
struct ProteinData {
  public:
    int gid;                 ///< unique global ID
    ProteinType property;    ///< intrinsic properties
    ProteinBindStatus bind;  ///< time-varying bind status
    double forceBind[2][3];  ///< spring force on bind MTs
    double torqueBind[2][3]; ///< spring torque on bind MTs

    /**
     * @brief Set default values for protein
     *
     * @return double
     */
    void setMockProtein() {
        bind.clear();
        property.setMockProteinType();
        std::fill(bind.pos, bind.pos + 3, 0);
    }

    /**********************************
     *
     * Get() const Functions
     *
     ***********************************/

    /**
     * @brief Get the pointer to doubel pos[3]
     *
     * @return const double*
     */
    const double *getPosPtr() const { return bind.pos; }

    /**
     * @brief Get gid
     *
     * @return int
     */
    int getGid() const { return gid; }

    /**
     * @brief Get BindingFactor from Unbound to Singly bound
     *
     * @param e
     * @param dt
     * @return double
     */
    double getBindingFactorUS(int e, double dt) const {
        assert(e == 0 || e == 1);
        assert(dt > 0);
        return property.ko_s[e] * property.Ka[e] * property.eps * dt;
    }

    /**
     * @brief Get UnbindingFactor from Singly bound to Unbound
     * TODO: Test in KMC_Testing
     *
     * @param e
     * @param dt
     * @return double
     */
    double getUnbindingFactorSU(int e, double dt) const {
        return property.ko_s[e] * dt;
    }

    /**
     * @brief Get BindingFactor from Singly to Doubly bound
     * TODO: Test in KMC_Testing
     *
     * @param e
     * @param dt
     * @return double
     */
    double getBindingFactorSD(int e, double dt) const {
        return property.ko_d[e] * property.Ke[e] * property.eps * dt;
    }

    /**
     * @brief Get UnbindingFactor from Doubly to Singly bound
     * TODO: Test in KMC_Testing
     *
     * @param e
     * @param dt
     * @return double
     */
    double getUnbindingFactorDS(int e, double dt, double KBT) const {
        // if (property.lambda == 0)
        return property.ko_d[e] * dt;
        // else {
        // double M = property.lambda * .5 * property.kappa *
        // SQR(getProteinForceLength() - property.freeLength) / KBT;
        // return property.ko_d[e] * dt * exp(M);
        //}
    }

    /**
     * @brief Get cut off calculation range when calculating 0 to 1 head bound
     *
     * @return double
     */
    double getRcutUS() const { return property.rc; }

    /**
     * @brief Get cut off calculation range when calculating 1 to 2 head bound
     * TODO: Test in KMC_Testing
     * @return double
     */
    double getRcutSD() const {
        // consistent tubuleDiameter as the original LUT construction
        double tubuleDiameter = property.LUTablePtr->getRodDiameter();
        // use the cutoff length as LUT construction
        return property.LUTablePtr->getNonDsbound() * tubuleDiameter;
    }

    /*! \brief Return the ID of the rod protein end is bound to
     *
     *
     * \param end Which protein end
     * \return ID of rod protein is bound to
     */
    int getBindID(int end) const { return bind.idBind[end]; }

    /**
     * @brief Get the ProteinForceLength, subtracting tubule diameter
     *
     * @return double the calculated length
     */
    /*
     *double getProteinForceLength() const {
     *    if (bind.idBind[0] == ID_UB || bind.idBind[1] == ID_UB) {
     *        return 0;
     *    } else {
     *        // consistent tubuleDiameter as the original LUT construction
     *        double tubuleDiameter = property.LUTablePtr->getRodDiameter();
     *        Evec3 r = ECmap3(bind.posEndBind[0]) - ECmap3(bind.posEndBind[1]);
     *        return (r.norm() - tubuleDiameter);
     *    }
     *}
     */

    /**
     * @brief Get the ProteinEndEndLength
     *
     * @return double the calculated length
     */
    /*
     *double getProteinEndEndLength() const {
     *    if (bind.idBind[0] == ID_UB || bind.idBind[1] == ID_UB) {
     *        return 0;
     *    } else {
     *        Evec3 r = ECmap3(bind.posEndBind[0]) - ECmap3(bind.posEndBind[1]);
     *        return r.norm();
     *    }
     *}
     */

    /**
     * @brief Get the reference to Protein Property
     *
     * @return const ProteinType&
     */
    const ProteinType &getProperty() const { return property; }

    /**
     * @brief get the ptr to LookupTable object
     *
     * @return LookupTable*
     */
    const LookupTable *getLUTablePtr() const { return property.LUTablePtr; }

    /**
     * @brief Get if the protein is walking along MT or not
     *
     * @return true
     * @return false
     */
    bool getWalkOrNot() const {
        return (bind.idBind[0] != ID_UB || bind.idBind[1] != ID_UB);
    }

    /**********************************
     *
     * Set() Functions
     *
     ***********************************/

    /**
     * @brief Set new position
     * Binding location is also jumped
     * @param newPos
     */
    /*
     *    void setPos(const double newPos[3]) {
     *        // current rvec
     *        Evec3 jump = ECmap3(newPos) - Emap3(bind.pos);
     *
     *        // set new pos
     *        std::copy(newPos, newPos + 3, bind.pos);
     *
     *        // move bind center if bind
     *        for (int e = 0; e < 2; e++) {
     *            if (bind.idBind[e] != ID_UB) {
     *                Emap3(bind.centerBind[e]) += jump;
     *            }
     *        }
     *        // update posEndBind and posProtein
     *        updateGeometryWithBind();
     *    }
     */

    /**********************************
     *
     * Update() Functions
     *
     ***********************************/

    /**
     * @brief update position by walking (S or D)
     *
     * @param dt
     */
    void updatePosWalk(double dt, double U01) {
        assert(getWalkOrNot());
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            // double v0 = calcEndWalkVelocity(0);
            // double v1 = calcEndWalkVelocity(1);
            double v0 = 0;
            double v1 = 0;
            bind.distBind[0] += v0 * dt;
            bind.distBind[1] += v1 * dt;
        } else if (bind.idBind[0] != ID_UB && bind.idBind[1] == ID_UB) {
            double v0 = property.vmax[0];
            bind.distBind[0] += v0 * dt;
        } else if (bind.idBind[0] == ID_UB && bind.idBind[1] != ID_UB) {
            double v1 = property.vmax[1];
            bind.distBind[1] += v1 * dt;
        } else {
            printf("walk while unbound is error\n");
            exit(1);
        }
        if (!property.walkOff) {
            updateWalkClampEnd();
        } else {
            updateWalkOffEnd(U01);
        }
        updateGeometryWithBind();
    }

    /**
     * @brief update bind when walk off
     *
     */
    void updateWalkOffEnd(double U01) {

        bool endOff[2] = {false, false};
        for (int end = 0; end < 2; end++) {
            if (bind.distBind[end] > bind.lenBind[end] * 0.5 ||
                bind.distBind[end] < bind.lenBind[end] * (-0.5)) {
                endOff[end] = true;
            }
        }

        if (endOff[0] && endOff[1]) {
            // both off, choose to unbind 1
            if (U01 < 0.5) {
                bind.setUnBind(0);
                bind.updatePosEndClamp(1);
            } else {
                bind.setUnBind(1);
                bind.updatePosEndClamp(0);
            }
        } else if (endOff[0]) {
            bind.setUnBind(0);
        } else if (endOff[1]) {
            bind.setUnBind(1);
        }
    }

    /**
     * @brief clamp distBind to end
     *
     */
    void updateWalkClampEnd() {
        bind.updatePosEndClamp(0);
        bind.updatePosEndClamp(1);
    }

    /**
     * @brief update position by diffusion (U)
     *
     * @param dt timestep size
     * @param N01x rngN01
     * @param N01y rngN01
     * @param N01z rngN01
     */
    void updatePosDiffuse(double dt, double N01x, double N01y, double N01z) {
        assert(!getWalkOrNot());
        const double diff = sqrt(dt * property.diffUnbound * 2);
        bind.pos[0] += diff * N01x;
        bind.pos[1] += diff * N01y;
        bind.pos[2] += diff * N01z;
    }

    /**
     * @brief update bind.posEndBind and bind.pos
     *
     */
    void updateGeometryWithBind() {
        bind.updatePosEndBind(0);
        bind.updatePosEndBind(1);
        bind.updatePosWithEndBind();
    }

    /**
     * @brief compute the spring binding force
     *
     * @return double
     */
    /*
     *void updateForceTorqueBind() {
     *    if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
     *        double tubuleDiameter = property.LUTablePtr->getRodDiameter();
     *        Evec3 r = Emap3(bind.posEndBind[0]) - Emap3(bind.posEndBind[1]);
     *        double length = r.norm();
     *        double force = (length - property.freeLength - tubuleDiameter) *
     *                       property.kappa;
     *        Evec3 f0 = -force * r.normalized();
     *        Evec3 f1 = -f0;
     *        // torque = r x f
     *        Evec3 torque0 =
     *            bind.distBind[0] * (Emap3(bind.directionBind[0]).cross(f0));
     *        Evec3 torque1 =
     *            bind.distBind[1] * (Emap3(bind.directionBind[1]).cross(f1));
     *        for (int k = 0; k < 3; k++) {
     *            forceBind[0][k] = f0[k];
     *            forceBind[1][k] = f1[k];
     *            torqueBind[0][k] = torque0[k];
     *            torqueBind[1][k] = torque1[k];
     *        }
     *    } else {
     *        // spring bind force is zero
     *        forceBind[0][0] = forceBind[0][1] = forceBind[0][2] = 0;
     *        forceBind[1][0] = forceBind[1][1] = forceBind[1][2] = 0;
     *        torqueBind[0][0] = torqueBind[0][1] = torqueBind[0][2] = 0;
     *        torqueBind[1][0] = torqueBind[1][1] = torqueBind[1][2] = 0;
     *    }
     *}
     */

    /**********************************
     *
     * Calc() const {} functions
     *
     ***********************************/

    /**
     * @brief compute the walking velocity of one end
     *
     * positive vmax towards plus end
     *
     * @param end 0 or 1
     * @return double computed walking velocity
     */
    /*
     *    double calcEndWalkVelocity(int end) const {
     *        assert(bind.idBind[end] != ID_UB);
     *
     *        double fproj =
     *            ECmap3(forceBind[end]).dot(ECmap3(bind.directionBind[end])) /
     *            abs(property.fstall);
     *        double vfrac = std::max(0.0, std::min(1.0, 1.0 + fproj));
     *        return property.vmax[end] * vfrac;
     *    }
     */

    /**********************************
     *
     * FDPS interface Functions
     *
     ***********************************/

    /**
     * @brief Get the Pos data field
     *  interface requried for FDPS
     * @return PS::F64vec3
     */
    // PS::F64vec3 getPos() const {
    //    return PS::F64vec(bind.pos[0], bind.pos[1], bind.pos[2]);
    //}
    const double *getPos() const { return bind.pos; }

    /**
     * @brief get neighbor search radius interface requried for FDPS
     *
     * @return double
     */
    double getRSearch() const { return property.freeLength * 10; }

    /**
     * @brief copyFromFP()
     *  interface requried for FDPS
     * @param fp
     */
    void copyFromFP(const ProteinData &fp) { *this = fp; }

    /**
     * @brief Set the Pos data field
     *  interface requried for FDPS
     * @param newPos
     */
    // void setPos(const PS::F64vec3 &newPos) {
    //    double newPos_[3] = {newPos.x, newPos.y, newPos.z};
    //    setPos(newPos_);
    //}
    // void setPos(const double *newPos) {
    //    double newPos_[3] = {newPos.x, newPos.y, newPos.z};
    //    setPos(newPos_);
    //}

    /**
     * @brief writeAscii file
     *  interface requried for FDPS IO
     * both cases can be read with the same routine
     * @param fptr
     */
    void writeAscii(FILE *fptr) const {
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            // protein has finite length
            // this should NOT out put nan
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.posEndBind[0][0], bind.posEndBind[0][1],
                    bind.posEndBind[0][2], //
                    bind.posEndBind[1][0], bind.posEndBind[1][1],
                    bind.posEndBind[1][2], //
                    bind.idBind[0], bind.idBind[1]);
        } else {
            // protein has zero length
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.idBind[0], bind.idBind[1]);
        }
    }

    /**********************************
     *
     * VTK format IO Functions
     *
     ***********************************/

    /**
     * @brief write VTK XML Parallel VTP file
     *
     * @tparam Container
     * @param protein
     * @param prefix
     * @param postfix
     * @param rank
     */
    //    template <class Container>
    //    static void writeVTP(const Container &protein, const int
    //    proteinNumber,
    //                         const std::string &prefix, const std::string
    //                         &postfix, int rank) {
    //        // each protein is a vtk polyline with two vertices and one cell

    //        // write VTP for basic data
    //        //  use float to save some space
    //        // point and point data
    //        // position always in Float64
    //        std::vector<double> pos(6 * proteinNumber);
    //        // label is for end 0 and end 1
    //        std::vector<float> label(2 * proteinNumber);

    //        // point connectivity of line
    //        std::vector<int32_t> connectivity(2 * proteinNumber);
    //        std::vector<int32_t> offset(proteinNumber);

    //        // protein data
    //        std::vector<int32_t> gid(proteinNumber);
    //        std::vector<int32_t> tag(proteinNumber);
    //        std::vector<int32_t> idBind(2 * proteinNumber);

    //#pragma omp parallel for
    //        for (int i = 0; i < proteinNumber; i++) {
    //            const auto &p = protein[i];
    //            // point and point data
    //            if (p.bind.idBind[0] != ID_UB && p.bind.idBind[1] != ID_UB) {
    //                // if both bind, posEndBind is valid
    //                pos[6 * i + 0] = p.bind.posEndBind[0][0];
    //                pos[6 * i + 1] = p.bind.posEndBind[0][1];
    //                pos[6 * i + 2] = p.bind.posEndBind[0][2];
    //                pos[6 * i + 3] = p.bind.posEndBind[1][0];
    //                pos[6 * i + 4] = p.bind.posEndBind[1][1];
    //                pos[6 * i + 5] = p.bind.posEndBind[1][2];
    //            } else {
    //                // else, shrink a point
    //                pos[6 * i + 0] = p.bind.pos[0];
    //                pos[6 * i + 1] = p.bind.pos[1];
    //                pos[6 * i + 2] = p.bind.pos[2];
    //                pos[6 * i + 3] = p.bind.pos[0];
    //                pos[6 * i + 4] = p.bind.pos[1];
    //                pos[6 * i + 5] = p.bind.pos[2];
    //            }

    //            label[2 * i] = 0;
    //            label[2 * i + 1] = 1;

    //            // connectivity
    //            connectivity[2 * i] = 2 * i;         // index of point 0 in
    //            line connectivity[2 * i + 1] = 2 * i + 1; // index of point 1
    //            in line
    //            // offset is the end of each line. in fortran indexing
    //            offset[i] = 2 * i + 2;

    //            // protein data
    //            // point data
    //            idBind[2 * i] = p.bind.idBind[0];
    //            idBind[2 * i + 1] = p.bind.idBind[1];
    //            // cell data
    //            gid[i] = p.gid;
    //            tag[i] = p.property.tag;
    //        }

    //        std::ofstream file(prefix + std::string("Protein_") + "r" +
    //                               std::to_string(rank) + std::string("_") +
    //                               postfix + std::string(".vtp"),
    //                           std::ios::out);

    //        IOHelper::writeHeadVTP(file);
    //        std::string contentB64; // data converted to base64 format
    //        contentB64.reserve(2 * proteinNumber);

    //        file << "<Piece NumberOfPoints=\"" << proteinNumber * 2
    //             << "\" NumberOfLines=\"" << proteinNumber << "\">\n";
    //        // Points
    //        file << "<Points>\n";
    //        IOHelper::writeDataArrayBase64(pos, "position", 3, file);
    //        file << "</Points>\n";
    //        // cell definition
    //        file << "<Lines>\n";
    //        IOHelper::writeDataArrayBase64(connectivity, "connectivity", 1,
    //        file); IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
    //        file << "</Lines>\n";
    //        // point data
    //        file << "<PointData Scalars=\"scalars\">\n";
    //        IOHelper::writeDataArrayBase64(label, "endLabel", 1, file);
    //        IOHelper::writeDataArrayBase64(idBind, "idBind", 1, file);
    //        file << "</PointData>\n";
    //        // cell data
    //        file << "<CellData Scalars=\"scalars\">\n";
    //        IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
    //        IOHelper::writeDataArrayBase64(tag, "tag", 1, file);
    //        file << "</CellData>\n";
    //        file << "</Piece>\n";

    //        IOHelper::writeTailVTP(file);
    //        file.close();
    //    }

    /**
     * @brief write VTK XML Parallel Header
     *
     * @param prefix
     * @param postfix
     * @param nProcs
     */
    //    static void writePVTP(const std::string &prefix, const std::string
    //    &postfix,
    //                          const int nProcs) {
    //        std::vector<std::string> pieceNames;

    //        std::vector<IOHelper::FieldVTU> pointDataFields;
    //        pointDataFields.emplace_back(1, IOHelper::IOTYPE::Float32,
    //        "endLabel"); pointDataFields.emplace_back(1,
    //        IOHelper::IOTYPE::Int32, "idBind");

    //        std::vector<IOHelper::FieldVTU> cellDataFields;
    //        cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
    //        cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "tag");

    //        for (int i = 0; i < nProcs; i++) {
    //            pieceNames.emplace_back(std::string("Protein_") +
    //            std::string("r") +
    //                                    std::to_string(i) + "_" + postfix +
    //                                    ".vtp");
    //        }

    //        IOHelper::writePVTPFile(prefix + "Protein_" + postfix + ".pvtp",
    //                                pointDataFields, cellDataFields,
    //                                pieceNames);
    //    }
};

static_assert(std::is_trivially_copyable<ProteinData>::value, "");
static_assert(std::is_default_constructible<ProteinData>::value, "");

/**
 * @brief FDPS writeAscii file header
 */
// class ProteinAsciiHeader {
//  public:
//    int nparticle;
//    double time;
//    void writeAscii(FILE *fp) const {
//        fprintf(fp, "%d \n %lf\n", nparticle, time);
//    }
//};

#endif
