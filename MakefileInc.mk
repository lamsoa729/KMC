# debug or not
DEBUGMODE:= yes

CXX= mpicxx
LINK= $(CXX)

# optimized. C++14 required by sctl
CXXFLAGS = -std=c++14 -fopenmp
LINKFLAGS= $(CXXFLAGS) 

RELEASE_OPT = -O3 -DNDEBUG
DEBUG_OPT = -O0 -g -DDEBUG

# EXTERNAL LIBRARIES one by one
# ALL STATIC LINK except mkl, openmp, and system runtime dynamic libs
# this can be different for each library
TRILINOSBASE = /usr/local
SFTPATH = /usr/local
# TRNG
USERINCLUDE += -I$(SFTPATH)/include/trng
USERLINK += $(SFTPATH)/lib/libtrng4.a
# EIGEN, header only
USERINCLUDE += -I$(SFTPATH)/include/eigen3
# Boost, header only
USERINCLUDE += -I/usr/local/Cellar/boost/1.69.0/include
# YAML
USERINCLUDE += -I$(SFTPATH)/include/
USERLINK += $(SFTPATH)/lib/libyaml-cpp.a

