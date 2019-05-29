#include "kmc_step_funcs.hpp"

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
    const double rPos = rod.getPos();
    // PS::F64vec3 rVec = rod.getPos();
    // double rPos[3] = {rVec[0], rVec[1], rVec[2]};
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
        pBind.setUnBind(head_bound);
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
        // PS::F64vec3 rVec = rod.getPos();
        const double rPos[3] = rod.getPos();
        // double rPos[3] = {rVec[0], rVec[1], rVec[2]};
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
