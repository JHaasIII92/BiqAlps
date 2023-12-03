#=============================================================================#
# This file is part of the Abstract Library for Parallel Search (ALPS).       #
#                                                                             #
# ALPS is distributed under the Eclipse Public License as part of the         #
# COIN-OR repository (http://www.coin-or.org).                                #
#                                                                             #
# Authors:                                                                    #
#                                                                             #
#          Yan Xu, Lehigh University                                          #
#          Aykut Bulut, Lehigh University                                     #
#          Ted Ralphs, Lehigh University                                      #
#                                                                             #
# Conceptual Design:                                                          #
#                                                                             #
#          Yan Xu, Lehigh University                                          #
#          Ted Ralphs, Lehigh University                                      #
#          Laszlo Ladanyi, IBM T.J. Watson Research Center                    #
#          Matthew Saltzman, Clemson University                               #
#                                                                             #
#                                                                             #
# Copyright (C) 2001-2019, Lehigh University, Yan Xu, Aykut Bulut, and        #
#                          Ted Ralphs.                                        #
# All Rights Reserved.                                                        #
#=============================================================================#


# $Id: Makefile.in 726 2006-04-17 04:16:00Z andreasw $

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = Biq
LBFGSBDIR = Lbfgsb.3.0

#
OBJDIR = ./objects/
# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  $(OBJDIR)BiqMain.o     \
	    $(OBJDIR)BiqModel.o    \
	    $(OBJDIR)BiqNodeDesc.o \
	    $(OBJDIR)BiqParams.o   \
	    $(OBJDIR)BiqSolution.o \
	    $(OBJDIR)BiqTreeNode.o \
	    $(OBJDIR)BiqUtil.o     \
		$(OBJDIR)BiqModelData.o

LBFGSBOBJS = $(OBJDIR)lbfgsb.o  \
			 $(OBJDIR)linpack.o \
			 $(OBJDIR)blas.o    \
	         $(OBJDIR)timer.o   

LBFGSBSRCS := lbfgsb.f linpack.f blas.f timer.f

ALLOBJS = $(OBJS) $(LBFGSBOBJS)
# CHANGEME: Additional libraries
ADDLIBS = -lopenblas  -lgfortran

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN package.                                                         #
##########################################################################

COIN_HAS_PKGCONFIG = TRUE
COIN_CXX_IS_CL = @COIN_CXX_IS_CL_TRUE@TRUE

# C++ Compiler command
CXX = g++
FC = gfortran
# C++ Compiler options
CXXFLAGS = -O2 -DNDEBUG 

ifeq ($(debug), 1)
    CXXFLAGS += -g
endif

FFLAGS =   -O -Wall -fbounds-check -Wno-uninitialized
# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/var/coin-or/dist/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/var/coin-or/dist/lib/pkgconfig: pkg-config --cflags alps`
#INCL = -I${prefix}/include/coin-or
INCL += $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/var/coin-or/dist/lib/pkgconfig: pkg-config --libs alps --static`
#LIBS = -L${exec_prefix}/lib -lAlps -lCoinUtils
LIBS += $(ADDLIBS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)



$(EXE): $(ALLOBJS)
	bla=;\
    for file in $(ALLOBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  $(LIBS) $(ADDLIBS) 

$(OBJDIR)%.o: %.cpp
	$(info Building $@ from $<)
	$(CXX) $(CXXFLAGS) $(INCL) -c $< -o $@

$(OBJDIR)%.o: %.f
	$(info Building $@ from $<)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJS) $(LBFGSBOBJS)


	
