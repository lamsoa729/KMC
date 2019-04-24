include ./MakefileInc.mk

# use DEBUG or RELEASE optimization flags
ifeq ($(DEBUGMODE), yes) 
CXXFLAGS += $(DEBUG_OPT)
LINKFLAGS+= $(DEBUG_OPT) 
else
CXXFLAGS += $(RELEASE_OPT)
LINKFLAGS+= $(RELEASE_OPT) 
endif

# inherit flags, includes and libs from Trilinos and pvfmm
include $(TRILINOSBASE)/include/Makefile.export.Trilinos
LINKFLAGS += $(Trilinos_EXTRA_LD_FLAGS)

# internal includes
SIMTOOLBOX = $(CURDIR)/SimToolbox
USERINCLUDE += -I$(SIMTOOLBOX) -I$(CURDIR) -I./SRC

# FDPS MPI and OpenMP
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL

# set of libraries to correctly link to the target
INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(USERINCLUDE)
LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS)
LIBRARIES = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(USERLINK)

# System-specific settings
SHELL = /bin/sh
SIZE =	size

# Files
SRCS = \
SRC/TubuleSystem.cpp \
SRC/SystemConfig.cpp \
Protein/ProteinConfig.cpp \
KMC/kmc.cpp \
KMC/lookup_table.cpp \
SimToolbox/Collision/CollisionSolver.cpp \
SimToolbox/Collision/CPSolver.cpp \
SimToolbox/Trilinos/TpetraUtil.cpp \
SimToolbox/Sylinder/Sylinder.cpp \
SimToolbox/Sylinder/SylinderConfig.cpp \
SimToolbox/Sylinder/SylinderSystem.cpp \
SimToolbox/Util/Base64.cpp \

# Definitions
EXE      := WinSim.X
TEST     := TEST.X
OBJ      := $(SRCS:.cpp=.o)
MAIN_OBJ := SRC/main.o
TEST_OBJ := SRC/test_main.o

all: $(EXE)

doc:
	doxygen TS-documentation.cfg

test: $(TEST)
	./$(TEST)

# Dependency files
OBJ_DEP = $(OBJ:.o=.d)
MAINOBJ_DEP = $(MAIN_OBJ:.o=.d)
TESTOBJ_DEP = $(TEST_OBJ:.o=.d)

# pull in dependency info for *existing* .o files
-include $(OBJ_DEP)
-include $(MAINOBJ_DEP)
-include $(TESTOBJ_DEP)

# Link rule
$(EXE): $(OBJ) $(MAIN_OBJ)
	$(LINK) $(OBJ) $(MAIN_OBJ)  -o $(EXE) $(LINKFLAGS) $(LIBRARY_DIRS) $(LIBRARIES)
	$(SIZE) $(EXE)

$(TEST): $(OBJ) $(TEST_OBJ)
	$(LINK) $(OBJ) $(TEST_OBJ)  -o $(TEST) $(LINKFLAGS) $(LIBRARY_DIRS) $(LIBRARIES)

%.o: %.cpp
	$(CXX) -MD -MP $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o


# remove compilation products
clean: 
	rm -vf $(EXE) $(TEST)  \
	$(OBJ) $(MAIN_OBJ) $(TEST_OBJ) \
	$(OBJ_DEP) $(MAINOBJ_DEP) $(TESTOBJ_DEP)
