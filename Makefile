# Copyright (C) 2003, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the four CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  cpp_example.o \
	MyNLP.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
#CXXFLAGS = -O2 -DNDEBUG -fPIC -fno-common -fexceptions
CXXFLAGS = -g -O0 -fsanitize=address -fPIC -fno-common -fexceptions

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

prefix=/Users/Juraj/switchdrive/Institution/usi/PhD/projects/Horenko/eSPA-code/ipopt-git/myinstall
exec_prefix=${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/Users/Juraj/switchdrive/Institution/usi/PhD/projects/Horenko/eSPA-code/ipopt-git/myinstall/lib/pkgconfig:/opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/bin/pkgconfig pkg-config --cflags ipopt` $(ADDINCFLAGS)
#INCL = -I${prefix}/include/coin-or   -DIPOPTLIB_BUILD $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/Users/Juraj/switchdrive/Institution/usi/PhD/projects/Horenko/eSPA-code/ipopt-git/myinstall/lib/pkgconfig:/opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/bin/pkgconfig pkg-config --libs ipopt`
#LIBS = -L${exec_prefix}/lib -lipopt -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/ /Users/Juraj/switchdrive/Institution/usi/PhD/projects/Horenko/eSPA-code/pardiso_mac/libpardiso700-MACOS-X86-64.dylib -lgfortran -lgomp -lquadmath -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lm  -ldl

all: $(EXE)

.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<
