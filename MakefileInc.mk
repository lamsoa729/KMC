# debug or not
DEBUGMODE:= yes

CXX= mpicxx
LINK= $(CXX)

# optimized. C++14 required by sctl
CXXFLAGS = -std=c++14 -fopenmp
LINKFLAGS= $(CXXFLAGS) 

RELEASE_OPT = -O3 -DNDEBUG
DEBUG_OPT = -O0 -g -DDEBUG

SFTPATH = $(HOME)/local
# Boost, header only
USERINCLUDE += -I$(SFTPATH)/include/boost/
# YAML
USERINCLUDE += -I$(SFTPATH)/include/
USERLINK += $(SFTPATH)/lib/libyaml-cpp.a

