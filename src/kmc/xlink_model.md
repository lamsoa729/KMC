# Introduction
The protein property implemented in this package is the two-stage explicit KMC model.
The reference for this model is:  
**Blackwell, R. et al. Physical determinants of bipolar mitotic spindle assembly and stability in fission yeast. Science Advances 3, e1601603 (2017).**  
Details on the physical protein model can be found in **Section 1.3 of the Supplementary Material** of this paper.
Note that there are many typos, errors, and inconsistencies in that document. 
In this markdown file, the notations are cleaned up, some symbols are changed, and the expressions to compute the four rates/probabilities are listed here. 

This document explains the coding convention and program behaviors, corresponding to the physical model.

# Units
All calculations in the code are based on \f$\mu m, pN, s\f$ units.
- \f$k_BT=0.00411 pN\cdot \mu m\f$ at room temperature 300K.
- MT diameter \f$D=0.025 \mu m = 25 nm\f$.
- Water viscosity \f$\eta = 0.00089 pN\cdot s \cdot \mu m^{-2} = 0.00089 Pa \cdot s\f$.
- Stress unit \f$1 Pa = 1 pN\cdot \mu m^{-2}\f$.

# The struct `ProteinType`
This is the struct holding the intrinsic properties of a protein. Different `ProteinType`s are initialized from the configuration file `ProteinConfig.yaml`.
Currently, the `ProteinType` does not change over time during the simulation. 
This can be made changeable in the future if necessary.
The code performs no check on the configuration. 
It is the user's duty to provide meaningful settings for `ProteinType`.

If `fixedEnd0: true` is set, the `freeNumber: 0` must be set, otherwise the code may crash, because a free protein cannot have a fixed and bound end. 
However, if `fixedEnd0: false` is set, `fixedLocationPerMT` can be arbitrarily set, which only means the initial position of each motor protein.

Variables in `ProteinType` are set by `ProteinConfig.hpp/cpp` routines by reading a `yaml` file named `ProteinConfig.yaml`. 
The per head and per protein properties are directly set as POD types.
The pointer to the lookup table `LUTablePtr` is set and refreshed after the MPI particle exchange stage and before the KMC stage within each timestep, according to the `tag` of each `ProteinType`, to make sure valid LUTables are used. 
The `LookupTable`s are centrally allocated and managed in the `TubuleSimulator` class.
The `LUTablePtr`s are merely users of those LookupTables. 

## The list of protein properties.
- **Lower case letters are reserved for protein properties**
- **Upper case letters are reserved for MT properties**
- **In real code, all dimensions will be converted to \f$\mu m\f$, \f$pN\f$, and \f$s\f$** 
### Per head (end) property
| Symbol    | Name                       | Number       | Dimension                       | Kinesin-5                      | `ProteinType` member variable |
| --------- | -------------------------- | ------------ | ------------------------------- | ------------------------------ | ----------------------------- |
| \f$K_a\f$     | Association constant       | \f$10\sim 100\f$ | \f$(\mu M)^{-1}\f$ per binding site | \f$90.9\f$                         | `Ka[2]` for two heads         |
| \f$K_E\f$     | Association constant       | \f$O(1)\f$       | dimensionless                   | \f$0.246 = 400/1625\f$             | `Ke[2]` for two heads         |
| \f$k_{o}^s\f$ | Singly bound turnover rate | \f$O(0.1)\f$     | \f$s^{-1}\f$                        | \f$0.11\f$ for                     | `ko_s[2]` for two heads       |
| \f$k_{o}^d\f$ | Doubly bound turnover rate | \f$O(0.05)\f$    | \f$s^{-1}\f$                        | half of \f$s\f$                    | `ko_d[2]` for two heads       |
| \f$v_m\f$     | Max velocity               | \f$10\sim1000\f$ | \f$nm/s\f$                          | \f$-50\f$ on polar aligned MT pair | `vmax[2]` for two heads       |

### Per protein property
| Symbol     | Name                                    | Number          | Dimension               | Kinesin-5 | `ProteinType` member variable |
| ---------- | --------------------------------------- | --------------- | ----------------------- | --------- | ----------------------------- |
| \f$\epsilon\f$ | Effective binding site density along MT | \f$125\sim 1650\f$  | site \f$\cdot \mu m^{-1}\f$ | \f$1625\f$    | `eps`                         |
| \f$\kappa\f$   | Spring constant                         | \f$O(1)\f$          | \f$pN/\mu m\f$              | \f$0.3\f$     | `kappa`                       |
| \f$\lambda\f$  | Unbinding load sensitivity              | \f$O(0 \sim 0.4)\f$ | dimensionless           | 0.25822   | `lambda`                      |
| \f$\ell_0\f$   | Rest length                             | \f$O(50)\f$         | \f$nm\f$                    | \f$53\f$      | `freeLength`                  |
| \f$r_c\f$      | Capture radius                          | **Note**        | \f$\mu m\f$                 | **Note**  | `rc`                          |
| \f$f_s\f$      | Stall force                             | \f$O(1 \sim 5)\f$   | \f$pN\f$                    | \f$5\f$       | `fstall`                      |
| \f$d_0\f$      | Diffusivity unbound                     | \f$O(1)\f$          | \f$\mu m^2 /s\f$            | \f$4.5\f$     | `diffUnbound`                 |

**Note**: The singly bound stage capture radius \f$r_c\f$ can be set in two different ways in ProteinConfig.yaml;
- \f$r_c = D/2 = 0.0125\f$. This makes the singly bound stage behavior the same as the old model.
- \f$r_c = (D + \ell_0)/2\f$. This makes the singly bound stage behavior follow the same convention as the new energy exponential factor including \f$-D\f$ in the doubly bound stage.


### Per protein behavior `bool` flags. 
- (To be implemented) Flag for binding anti-parallel MT only `bool bindAntiParallel`.
- (To be implemented) Flag for one end fixed `bool fixedEnd0`.
- Flag for walking off the ends `bool walkOff`.


# The struct `ProteinBindStatus`
This is the struct holding the position and binding (dynamic) information of a protein. 
Each protein is tracked as a single point if at unbound or singly bound state (U or S), but tracked as a spring with two ends if at doubly bound state (D).

**Important** This is the only data that is updated by the KMC calculations.

Therefore, the datafields `double pos[3]` and `double posEndBind[2][3]` should be carefully handled.
More on this in the `ProteinData` section.

The definitions of default values are:
```cpp
constexpr int ID_UB = -1;
constexpr double NAND = std::numeric_limits<double>::quiet_NaN();
```
Here is a list of `ProteinBindStatus` data fields.  
| variable                     | default value                         | description                                                                                                  |
| ---------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `double pos[3]`              | `{NAND,NAND,NAND}`                    | The (lab frame) center position of this protein. Always valid                                                |
| `bool changeBind[2]`         | `{true, true}`                        | If `false`, update binding information only and bypass the KMC step.                                         |
| `int idBind[2]`              | `{ID_UB, ID_UB}`                      | The gid of each bound MT. Set to `ID_UB` if unbound                                                          |
| `int rankBind[2]`            | `{ID_UB, ID_UB}`                      | The MPI rank of each bound MT. \f$\in[0,nProcs-1]\f$, otherwise set to `ID_UB`.                                  |
| `double lenBind[2]`          | `{NAND, NAND}`                        | The length of each bound MT.                                                                                 |
| `double distBind[2]`         | `{NAND, NAND}`                        | The distance to the center of each bound MT. \f$\in [-lenBind/2,lenBind/2]\f$, positive is towards the plus end. |
| `double posEndBind[2][3]`    | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (lab frame) location of each protein end (head) when bound.                                              |
| `double centerBind[2][3]`    | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (lab frame) location of the center of each bound MT                                                      |
| `double directionBind[2][3]` | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (normalized) orientation vector of each bound MT                                                         |

# The Binding/Unbinding probabilities and rates in KMC
The KMC rates and probabilities depend on both properties and geometry.
In the following, the superscript \f$i=0,1\f$ denotes the two heads of a protein. 
For the geometric dependent part, uppercase symbols are for MT, like \f$L,D\f$ for length and diameter.
Lowercase symbols are for protein, like \f$\ell,r_c\f$ for length and capture radius.

## Between Unbound and Singly Bound Stages
| Probability                              | Rate                                                                                                                                |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| \f$P_{on}^{s,i} = k_{on}^{s,i}\Delta t\f$    | \f$k_{on}^{s,i} = \left( \frac{L_c}{2r_c} \right) \left(\frac{3 K_a^i }{4\pi r_c^3} \right) \left( 2r_c\epsilon\right) k_{o}^{s,i}\f$ |
| \f$P_{off}^{s,i} = k_{off}^{s,i} \Delta t\f$ | \f$k_{off}^{s,i}=k_{o}^{s,i}\f$                                                                                                         |
\f$L_c\f$ is the length of MT within the capture-sphere with radius \f$r_c\f$ of protein.
\f$\Delta t\f$ is the discretized timestep size.

## Between Singly Bound and Doubly Bound Stages
| Probability                                            | Rate                                                                                                                                                                              |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| \f$\displaystyle P_{on}^{d,i} = k_{on}^{d,i} \Delta t\f$   | \f$\displaystyle k_{on}^{d,i} =  \left\{D\epsilon K_E^i \int\exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(\ell(s)-\ell_0-D\right)^2\right] \frac{1}{D} ds \right\}  k_{o}^{d,i}\f$ |
| \f$\displaystyle P_{off}^{d,i} = k_{off}^{d,i} \Delta t\f$ | \f$\displaystyle k_{off}^{d,i} = \exp\left[\frac{\lambda \kappa \beta}{2} \left(\ell-\ell_0-D\right)^2\right] k_{o}^{d,i}\f$                                                          |

## Protein Force vs Protein Length
In the model, each protein is tracked by either a point (U or S stage) or a spring (D stage).
In the U or S stage, there is no confusion because the protein is a point and force is always zero.
In the D stage, however, we make a distinction between the 'geometric length' and the 'force length' for a protein.

The geometric length is simply the length \f$\ell\f$ between two `posEndBind[3]` along both MT **center line**s.
The 'force length' is defined as \f$\ell - D\f$, where \f$D\f$ is the MT diameter.
Here the subtraction off the MT diameter \f$D\f$ is because the length \f$\ell\f$ is computed between centerlines, but the actual binding happens on the MT surfaces. 
Subtracting off \f$D\f$ is an easy and crude approximation for this effect.
This also significantly softens the 'contractile' effects of proteins.

The stretch of a protein is \f$\ell-\ell_0-D\f$, the force of the protein is \f$\kappa (\ell-\ell_0 - D)\f$, and the energy is \f$\frac{1}{2}\kappa (\ell-\ell_0-D)^2\f$.

**IMPORTANT**
Under this model, a protein can also generate extensile stress if it is compressed.

## The (dimensionless) S\f$\rightarrow\f$D probability integral
\f[\int_{s_A}^{s_B}\exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(\ell(s)-\ell_0-D\right)^2\right] \frac{1}{D} ds\f]

The bound head \f$\alpha\f$ of the protein is already fixed at a point somewhere in space (another MT). 
We first compute the perpendicular foot from head \f$\alpha\f$ to the infinite centerline, where the distance is denoted as \f$\ell_{m}\f$, meaning the 'minimal distance' from \f$\alpha\f$ to the infinite centerline. 
\f$\ell_m\f$ corresponds to the variable `distPerp` in the code.

We then put an axis towards the plus direction of the MT on the MT centerline, where \f$s=0\f$ marks the perpendicular foot corresponding to \f$\ell_m\f$.
The tubule is marked by bounds \f$s_A\f$ and \f$s_B\f$ on its infinite centerline. 
\f$s_B-s_A=L\f$, the MT length. 
Depending on the geometry, there are in total three cases:
- \f$0<s_A<s_B\f$
- \f$s_A<0<s_B\f$
- \f$s_A<s_B<0\f$

The integral is then non-dimensionalized with MT diameter \f$D\f$:
- \f$s^\prime = s /D\f$
- \f$\ell^\prime = \ell /D\f$

The reason why choosing \f$D\f$ as the length scale is because
- Different proteins with different \f$\ell_0\f$ exists in the system. Choosing \f$D\f$ gives a consistent basis for all proteins.
- Actual stretching of proteins happens on the length scale of \f$D\f$ (collision, alignment, etc.), instead of \f$\ell_0\f$.

The integral becomes
\f[\int_{s_{A}^\prime}^{s_{D}^\prime}\exp\left[-M\left(\sqrt{{\ell_m^\prime}^2 +{s^\prime}^2}-\ell_0^\prime-1\right)^2\right] ds^\prime\f]
A dimensionless parameter \f$M\f$ can be identified, as the binding free energy factor:
\f[M=\frac{(1-\lambda)\kappa\beta D^2}{2}\f]
In general, \f$M\in(0.1,10)\f$ for proteins used in this model.

This nondimensionalization avoids numerical integration of a sharp peak over length scale \f$L\f$. 
Instead, the integration of length scale \f$D\f$ is soft and smooth, straightforward to be integrated by Gauss-Kronrod or other 1D quadrature schemes. 

The behavior of the integrand depends on the minimal distance \f$\ell_m^\prime\f$:
- \f$\ell_m^\prime > 1+\ell_0^\prime\f$: the integrand is single peaked at \f$s^\prime=0\f$. The proein is always stretched.
- \f$\ell_m^\prime < 1+\ell_0^\prime\f$: the integrand is double peaked at \f$s^\prime=\pm\sqrt{(\ell_0^\prime+1)^2-\ell_m^{\prime 2}}\f$. The protein can be compressed or stretched depending on \f$s\f$.

Due to symmetry about \f$s^\prime=0\f$, the integral is tabulated as follows:
\f[\int_0^{s_{bound}^\prime}\exp\left[-M\left(\sqrt{{\ell_m^\prime}^2 +{s^\prime}^2}-\ell_0^\prime-1\right)^2\right] ds^\prime\f]
Where \f$M\f$ and \f$\ell_0^\prime\f$ are fixed for each species. The integral is tabulated over \f$\ell_m^\prime\f$ and \f$s_{bound}^\prime\f$ as a 2D regular table. The tabulation is truncated at the integrand \f$\exp\left[-M\left(\sqrt{{\ell_m^\prime}^2 +{s^\prime}^2}-\ell_0^\prime-1\right)^2\right]=eps\f$. This gives the bounds:
- \f$0<s_{bound}^\prime<\sqrt{l_{UB}^{\prime 2}-1}\f$.
- \f$1<\ell_m^\prime<l_{UB}^\prime\f$.

Here \f$\ell_{UB}^\prime = \ell_0^\prime +1 +\sqrt{-\frac{1}{M}\ln eps}\f$, is the (dimensionless) truncation length corresponding to \f$eps\f$. 
The table can be constructed as a \f$32\times64\f$ matrix. \f$32\f$ grids along \f$\ell_m^\prime\f$ and \f$64\f$ grids along \f$s_{bound}^\prime\f$. 
Lookup entries lower than the bounds generate an error, while values higher than the bounds return zero.


## Dimensionless groups for binding/unbinding process
Identifying the dimensionless groups helps understanding the effects of the numerical parameter values.
| Name                                                                                                                          | Kinesin-5 | Explanation                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------- | --------- | ------------------------------------------------------------------------------------------ |
| \f$\frac{3 K_a}{4\pi r_c^3}\f$                                                                                                   | ?         | Dimensionless association constant                                                         |
| \f$2\epsilon r_c\f$                                                                                                               | ?         | Total number of binding sites in capture radius                                            |
| \f$\frac{L_c}{2r_c}\f$                                                                                                           | \f$[0,1]\f$   | Geometric capturing factor depending on the configuration of MP and MTs                    |
| \f$\left(\frac{3 K_a}{4\pi r_c^3} \right) \left( 2 \epsilon r_c \right)\f$                                                       | ?         | Ratio between on and off rate                                                              |
| \f$\left( \frac{L_c}{2 r_c} \right) \left(\frac{3 K_a^i }{4\pi r_c^3} \right) \left( 2r_c\epsilon\right) k_{o}^{s,i}\Delta t\f$ | ?         | Scale of binding probability per timestep. Should be \f$<1\f$, otherwise \f$\Delta t\f$ too large. |
| \f$D\epsilon K_E^i\f$                                                                                                             | ?         | Dimensionless association constant                                                         |
| \f$M=\frac{(1-\lambda)\kappa\beta D^2}{2}\f$                                                                                      | ?         | Scale of doubly bound unbinding rate.                                                      |

# The struct `ProteinData`
## Definition
This struct describes one protein, for all different types.
This is the data structure used in actual dynamic simulations. 
Each protein has one unique non-negative integer global ID `gid`, a `ProteinProperty` struct, a `ProteinBindStatus` struct, a `forceBind[2][3]`, and a `torqueBind[2][3]`.

## Member functions
The member functions (besides simple `set` and `get` functions) are grouped into two categories:
- `updateXXX()` functions update the internal member variables. No return values (`void updateXXX()`).
- `calcXXX() const` functions calculate something with member variables and some given non-member information. Internel variables are kept constant, and the calculated values are returned.

## Initialization
If the initial protein configuration is not provided as an input file `ProteinInitial.txt`, initial configurations will be generated according to the configuration file `ProteinConfig.yaml`. 
Each protein will be sequentially numbered and generated on mpi rank 0 after the maximum `gid` of sylinders.
The position of each protein, `pos[3]`, is either generated from uniform random numbers inside the initBox, or set according to tubule configuration if `fixedEnd0: true` is set in `ProteinConfig.yaml`.

If the protein initial configuration file `ProteinInitial.txt` is provided, either from a history file of the previous simulation, or hand-written, the `gid` will be read from the configuration file. 
*No check* on this `gid` is performed. 
It is the user's job to provide meaningful data.

## Protein position in simulation
There are two scenarios where the position of a protein is read or set: (1) inside FDPS routines and (2) outside FDPS routines.  

Inside FDPS routines, the `ProteinData` struct appears in two cases: as a `FDPS::FullParticle` or a `FDPS::EssentialParticleI`. 
- When `ProteinData` is used as a `FDPS::FullParticle`, `void setPos(const PS::F64vec3 &)` is called **ONLY** inside the `adjustPositionIntoRootDomain()` function, to set the position according to the periodic boundary condition along a certain direction of a box.
- When `ProteinData` is used as a `FDPS::EssentialParticleI`, `void setPos(const PS::F64vec3 &)` is **NEVER** called. This conclusion also applies to the case when `ProteinData` is used as `EPI` in `MixedPairInteration`. In this case, `PS::F64vec3 getPos() const` is called, and must return a meaningful position in the (possibly periodic) simulation box `RootDomain`, otherwise `FDPS` may work abnormally and is hard to detect at runtime.

Therefore, when `void setPos(const PS::F64vec3 &newPos)` is called, the jump (`newPos-pos`) must represent periodic jumps across periodic images of the simulation box. 
The `double centerBind[2][3]`, if not `NAND`, should also be moved as the same jump `newPos-pos`. 
The implicit assumption is that there is always a MT across the periodic image jump there, so the `double centerBind[2][3]` remains always valid after the jump.
This is necessary to make FDPS work.

Outside FDPS routines, the location is accessed by overloaded functions `const double * getPosPtr() const` and `void setPos(const pos[3])`.
These functions get and set the `pos` of the protein, with the same behavior as the FDPS interface versions.

## Protein cycle within one timestep

At each time step, the following events happen **in the given order**.

| Protein                    | MT                   | Comment                                                              |
| -------------------------- | -------------------- | -------------------------------------------------------------------- |
|                            | `prepareStep()`      | partition, exchange, update map, etc                                 |
| apply PBC                  |                      | `FDPS::adjustPositionIntoRootDomain()`                               |
| mpi exchange               |                      | `FDPS::exchangeParticle()` accoridng to `domainInfo` of MT partition |
| run KMC step               |                      | `calcBindInteraction()`, which calls `MixedPairInteraction`.         |
| `updateGeometryWithBind()` |                      |                                                                      |
| walk or diffuse            |                      | if walk, check `walkOff`                                             |
| `updateForceTorqueBind()`  |                      |                                                                      |
| save protein info to file  |                      | save to XML vtk file                                                 |
| apply force & torque to MT | `setForceNonBrown()` | `ZDD` and `CommMPI` classes are involved                             |
|                            | `runStep()`          | MT history file is saved inside `runStep()` before MT movement       |

**Note**
- When `walkOff == true`, protein unbinds when it walks past the ends of bound MT. If two ends simultaneously walk past the minus ends, one end is chosen (with equal probability for both ends) to walk off and the other end is clamped to the end and is still bound. This protein then continues the normal KMC cycle in the next timestep 

<!-- ### Pre-KMC
The property `fixedEnd0` is handled here. If `true`, the `changeBind[0]=false` is set to bypass the KMC step for this end. 
In this case, the protein switches between singly and doubly bound states only. 
In the `CalcProteinBind()` step, the information of the fixed MT is updated because MTs are moving. 
### KMC (near-neighbor)
There are two cases: `changeBind[2]={false, true}` or `changeBind[2]={true,true}`. In the first case, `KMC_1_u` or `KMC_2_1` is executed depending on the binding status. In the second case, `KMC_0_1`, `KMC_1_u`, or `KMC_2_1` is executed depending on the binding status.
### Post-KMC
After this check, the protein motion is moved according to given diffusivity and walking velocity. 
Finally, the protein force is calculated and added to the tubules by using the `ZDD` and `CommMPI` classes. The protein binding stress is also calculated here. -->

# The struct `ProteinConfig`
This is a class reading in a `yaml` file to specify protein properties and initial configuration. 
An example config `yaml` file is included, named `ProteinConfig_example.yaml`.

The file specified all protein types as a `yaml` list:
```yaml
proteins:
  - ....# specify protein type 0
  - ....# specify protein type 1
  - ....# specify protein type 2
  - .... 
```
Currently there is no limit of how many types of proteins are specified.

Each type of protein includes a property section:
```yaml
    tag: 0 # Type 0,
    #properties:
    walkOff: true
    bindAntiParallel: false
    fixedEnd0: false
    freeLength: 0.08 # um
    kappa: 1.0 # pN/um
    fstall: 1.0 # pN
    lambda: 0 # dimensionless
    diffUnbound: 5.0 # um^2/s
    vmax: [1.0, 1.0] # um/s, positive towards plus end
    # KMC parameters
    eps: 1625 # 13 protofilaments at 8nm per block
    Ka: [53.05, 10.0] # (uM/L)^{-1}
    ko_s: [0.003916, 0.001958] # 1/s
    Ke: [0.246, 0.1] # dimensionless
    ko_d: [0.003916, 0.001958] # 1/s
```
and an initial configuration section:
```yaml
    freeNumber: 0
    fixedLocationPerMT: [-2,0,1] # given in [-1,1]. otherwise random
``` 
The order of properties and configurations specified in this file does not matter.
Most properties are self-explainatory. 
A single parameter refers to the property per protein, and an array of two parameters refers to the properties of each head (end). 
These are read and stored in the `ProteinProperty` struct and every protein has an independent copy of all the properties.

The special property `tag` is an arbitrary integer given by the user to mark the types. The tags do not have to be sequentially ordered.
For example, the tags can be specified like this:
```yaml
proteins:
  - tag: 5 # for Kinesin-5
    .....  # other Kinesin-5 properties
  - tag: 1 # for Kinesin-1
    .....  # other Kinesin-1 properties
``` 

The initial configuration section specifies specifies the free and bound proteins in the beginning.
`freeNumber` denotes the number of proteins that are uniformly and randomly distributed in the simulation box, regardless of whether periodic boundary conditions are set.
`fixedLocationPerMT` is an array and denotes the number and locations of proteins singly bound to each MT in the beginning. 
An entry in this array specifies the location.
A number \f$\in[-1,1]\f$ specifies the location, from minus end (-1) to positive end (1).
A number out of that range (e.g. 2, -2, 3, etc) specifies that a random location is chosen. 
For example:
- `fixedLocationPerMT: []` means no initially bound protein.
- `fixedLocationPerMT: [-3,-2,0,1]` means 4 proteins (of this type) are bound to each MT in the beginning.
One at the center (0), one at the positive end (1), two are randomly chosen (-3,-2). 

# Protein IO
## Data output
Output frequency is controlled by 
```yaml
dt: 0.001 # s
timeTotal: 0.1 # s
timeSnap: 0.001 # s
```
in `RunConfig.yaml` file. 
This example will generate 1 (`timeSnap` divide by `dt`) snapshots per timestep, in total 100 (`timeTotal` divide by `timeSnap`) snapshots.

Each snapshot includes one ASCII file (`.dat`) and a group of binary (base64) XML VTK file (`.pvtp` and `.vtp`). 
### ProteinAscii_NUMBER.dat
This is the ASCII file to record protein geomery and bound MT ID only.
`NUMBER` is sequentially ordered from 0 to mark the sequence of files.
All proteins from all mpi ranks are combined into a single file.
Each protein takes a line, and is exported as
```c
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
```
**NOTE** The order of protein `gid` appears in this dat file is **NOT** specified and is **NOT** fixed.
### Protein_NUMBER.pvtp and Protein_rX_NUMBER.vtp 
Each mpi rank writes its own partition of protein data. `rX` marks the mpi rank. The data can be read and visualized by `Paraview` using the included `TubuleVis.pvsm` file for `Paraview` higher than version 5.6.

The explanation of binary (base64) XML VTK file cna be found on VTK official website.

## Data input
The file `ProteinAscii_NUMBER.dat` can be copied to `ProteinInitial.dat` to be read by the main executable program to set the initial configuration and binding status of proteins. 
Binding information will be reconstructed with proper `TubuleInitial.dat`, `RunConfig.yaml`, and `ProteinConfig.yaml`.

**Note**
- The protein tags in `ProteinConfig.yaml` must include all tags appearing in `ProteinInitial.dat`. Properties of the same tag can be different from a previous run.
- The `TubuleInitial.dat` must match the `ProteinInitial.dat` if any protein is in S or D bound stage.  
- The simbox geometry in `RunConfig.yaml` must match if any protein is in S or D bound stage.  
