/**
 * @file ProteinType.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef PROTEINTYPE_HPP_
#define PROTEINTYPE_HPP_

#include "lookup_table.hpp"

#include <cstdio>

/**
 * @brief Specify protein behavior
 *
 */
struct ProteinType {
  public:
    // per protein properties
    int tag = 0;              ///< user-assigned integer tag for different types
    bool walkOff = true;      ///< walf off the end, i.e., no 'end-dewelling'
    double fixedLocation = 0; ///< in [-1,1]
    double freeLength = 0.05; ///< the free length. unit um
    double kappa = 204.7; ///< Spring constant when attached to MT. unit pN/um
    double fstall = 1.0;  ///< stall force. unit pN
    double diffUnbound = 4.5; ///< Unbounded diffusivity, unit um^2/s
    double lambda = 0;        ///< dimensionless unbinding load sensitivity

    // TODO: Implement these two functions
    bool bindAntiParallel = false; ///< only bind anti parallel MTs
    bool fixedEnd0 = false; ///< end0 is fixed to a MT at the given location

    // per head properties
    double vmax[2]; ///< singly bound max velocity for each end. unit um/s

    // kMC data field
    double rc;  ///< the capture radius of protein (user set)
    double eps; ///< Number of crosslinker binding sites per um of MT
    double ko_s[2];
    ///< Turnover rate for each end for KMC_u_s and KMC_s_u steps
    double ko_d[2];
    ///< Turnover rate for each end for KMC_s_d and KMC_d_s steps
    double Ka[2];
    ///< Association constant when neither head is  bound
    double Ke[2];
    ///< Equilibrium constant when other end is bound

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

    void setMockProteinType() {
        Ka[0] = Ka[1] = 1.0;
        Ke[0] = Ke[1] = 1.0;
        ko_s[0] = ko_s[1] = 1.0;
        ko_d[0] = ko_d[1] = 1.0;
        vmax[0] = vmax[1] = 1.0; // max velocity. unit um/s
        kappa = 1.0;             // Spring constant when attached to MT
        eps = 1.0;
        rc = .5;          // capture radius
        freeLength = 1.0; // the 'bind search' radius
        fstall = 1.0;     ///< stall force. unit pN
        tag = 0;          ///< user-assigned integer tag for different types
        walkOff = true;   ///< walf off the end, i.e., no 'end-dewelling'
        lambda = 0;
    }
};

#endif
