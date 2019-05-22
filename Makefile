include ./MakefileInc.mk

# use DEBUG or RELEASE optimization flags
ifeq ($(DEBUGMODE), yes) 
CXXFLAGS += $(DEBUG_OPT)
LINKFLAGS+= $(DEBUG_OPT) 
else
CXXFLAGS += $(RELEASE_OPT)
LINKFLAGS+= $(RELEASE_OPT) 
endif

# set of libraries to correctly link to the target
INCLUDE_DIRS = -I$(CURDIR) -I$(CURDIR)/KMC $(USERINCLUDE) 
LIBRARIES = $(USERLINK)

# System-specific settings
SHELL = /bin/sh
SIZE =	size

# Files
SRCS = \
	KMC/kmc.cpp \
	KMC/lookup_table.cpp \
	KMC/integrals.cpp

# Definitions
TEST     := TEST.X
OBJ      := $(SRCS:.cpp=.o)
TEST_OBJ := tests/test_main.o

all: $(EXE)

doc:
	doxygen kmc-documentation.cfg

test: $(TEST)
	./$(TEST)

# Dependency files
OBJ_DEP = $(OBJ:.o=.d)
TESTOBJ_DEP = $(TEST_OBJ:.o=.d)

# pull in dependency info for *existing* .o files
-include $(OBJ_DEP)
-include $(TESTOBJ_DEP)

# Link rule
$(TEST): $(OBJ) $(TEST_OBJ)
	$(LINK) $(OBJ) $(TEST_OBJ)  -o $(TEST) $(LINKFLAGS) 

%.o: %.cpp
	$(CXX) -MD -MP $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o


# remove compilation products
clean: 
	rm -vf $(EXE) $(TEST)  \
	$(OBJ) $(TEST_OBJ) \
	$(OBJ_DEP) $(TESTOBJ_DEP)
