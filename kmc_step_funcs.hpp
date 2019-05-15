#ifndef KMC_STEP_FUNCS_HPP
#define KMC_STEP_FUNCS_HPP

#include "helpers.hpp"
#include "kmc.hpp"
#include "lookup_table.hpp"
#include "macros.hpp"
#include <array>
#include <cassert>
#include <cmath>

#include "Protein/ProteinBindStatus.hpp"
#include "Protein/ProteinData.hpp"
#include "Protein/ProteinType.hpp"
#include "SimToolbox/Sylinder/SylinderNear.hpp"
//#include "SylinderNear.hpp"
//#include "example_objs/ExampleRod.hpp"
// template class KMC<ExampleRod>;

// These can be changed depending on program and desired behavior

// KMC step for unbound protein
void KMC_U(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double roll,
           ProteinBindStatus &pBind);

// KMC step for singly bound protein
void KMC_S(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
           double rollVec[3], ProteinBindStatus &pBind);

// KMC step for doubly protein
void KMC_D(const ProteinData &pData, const int Npj,
           const SylinderNearEP *const *ep_j,
           const std::vector<int> &uniqueFlagJ, double dt, double KBT,
           double roll, ProteinBindStatus &pBind);

#endif /* KMC_STEP_FUNCS_HPP */
