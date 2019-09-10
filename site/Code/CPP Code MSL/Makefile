# The generated "Makefile" from "Makefile.cmake" is only usable after
# the Armadillo library has been configured and installed by CMake.

TMW_ROOT=/n/sw/matlab-R2011a

KNDIR=/n/sw/knitro-6.0.0-z
KNRELEASE=600
KNLIB_STATIC = $(KNDIR)/lib/libknitro$(KNRELEASE).a
KNLIB_DYNAMIC = $(KNDIR)/lib/libknitro.so
KFLAGS=-I$(KNDIR)/include

MFLAGS= -I$(TMW_ROOT)/extern/include
MFLAGS2= -I$(TMW_ROOT)/extern/include -largeArrayDims
MLIBS= -L$(TMW_ROOT)/bin/glnxa64 -lmat -lmx
MEXFLAGS= -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread -DNDEBUG
MEXLIB = -pthread -shared -Wl,--version-script,$(TMW_ROOT)/extern/lib/glnxa64/mexFunction.map  -lmex

AFLAGS=-I/usr/local/include
ALIBS=-L/usr/local/bin -larmadillo

CXX=g++
#CXX=g++
## Under MacOS you may have an old compiler as default (e.g. GCC 4.0).
## However, GCC 4.2 or later is available and preferable due to better
## handling of template code.

#CXX=CC
## When using the Sun Studio compiler


# flags configured by CMake

#EXTRA_LIB_FLAGS = -library=sunperf
## When using the Sun Studio compiler


ifeq (true,true)
  BOOST_INCLUDE_FLAG = -I /opt/local/include
endif



LIB_FLAGS = -larmadillo $(EXTRA_LIB_FLAGS) $(MLIBS) $(KNLIB_DYNAMIC)
## NOTE: on Ubuntu and Debian based systems you may need to add 
## -lgfortran to LIB_FLAGS



OPT = -O3 -funroll-all-loops   -fopenmp -m64
## As the Armadillo library uses recursive templates,
## compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: produces programs which achieve most of the possible speedup
## -O3: produces programs which have almost all possible speedups,
##      but compilation takes considerably longer


#OPT = -xO4 -xannotate=no
## When using the Sun Studio compiler


# EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling 
## all source files into one program in a single hit.


#DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level
## debugging.  Lots of debugging information will
## be printed when a compiled program is run.
## Please enable this option when reporting bugs.


# FINAL = -DARMA_NO_DEBUG 
## Uncomment the above line to disable Armadillo's checks.
## DANGEROUS!  Not recommended unless your code has been
## thoroughly tested.


#
#
#

CXXFLAGS = $(BOOST_INCLUDE_FLAG) $(MFLAGS) $(KFLAGS) $(AFLAGS) $(FINAL) $(OPT) $(EXTRA_OPT) -fPIC 
CXX2FLAGS = $(BOOST_INCLUDE_FLAG) $(MFLAGS2) $(KFLAGS) $(AFLAGS) $(FINAL) $(OPT) $(EXTRA_OPT) -fPIC 

OBJS= armahelp.o globals.o msllcpp.o

all:  mslestimate choiceprob test llevalmex

mslestimate: mslestimate.o armahelp.o globals.o msllcpp.o
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJS)  $<  $(LIB_FLAGS)

mslconstr: mslconstr.o armahelp.o globals.o msllcpp.o
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJS)  $<  $(LIB_FLAGS)

test: test.o armahelp.o globals.o msllcpp.o
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJS)  $<  $(LIB_FLAGS)

choiceprob: choiceprob.o globals.o msllcpp.o armahelp.o
	$(CXX) -fopenmp -fPIC -o choiceprob.mexa64 $(OBJS)  $< -larmadillo $(EXTRA_LIB_FLAGS) $(MLIBS) $(MEXLIB)

llevalmex: llevalmex.o globals.o msllcpp.o armahelp.o
	$(CXX) -fopenmp -o lleval.mexa64 $(OBJS)  $< -larmadillo $(EXTRA_LIB_FLAGS) $(MLIBS) $(MEXLIB)

choiceprob.o: choiceprob.cpp globals.cpp   
	$(CXX) $(CXX2FLAGS) $(MEXFLAGS) -c choiceprob.cpp

llevalmex.o: llevalmex.cpp globals.cpp   
	$(CXX) $(CXX2FLAGS) $(MEXFLAGS) -c llevalmex.cpp

globals.o: globals.cpp
	$(CXX) $(CXXFLAGS)  -c globals.cpp

test.o: test.cpp globals.cpp
	$(CXX) $(CXXFLAGS)  -c test.cpp

mslestimate.o: mslestimate.cpp globals.cpp
	$(CXX) $(CXXFLAGS)  -c mslestimate.cpp

mslconstr.o: mslconstr.cpp globals.cpp
	$(CXX) $(CXXFLAGS)  -c mslconstr.cpp

msllcpp.o : msllcpp.cpp globals.cpp
	$(CXX) $(CXXFLAGS)  -c msllcpp.cpp

armahelp.o: armahelp.cpp
	$(CXX) $(CXXFLAGS)  -c armahelp.cpp	


.PHONY: clean

clean:
	rm -f test *.o lleval.mexa64 mslconstr


