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

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  BiqMain.o \
	BiqModel.o \
	BiqNodeDesc.o \
	BiqParams.o \
	BiqSolution.o \
	BiqTreeNode.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = /home/james/Alps/examples/Biq
VPATH = /home/james/Alps/examples/Biq

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN package.                                                         #
##########################################################################

COIN_HAS_PKGCONFIG = TRUE
COIN_CXX_IS_CL = @COIN_CXX_IS_CL_TRUE@TRUE

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O2 -DNDEBUG 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/james/dist/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/home/james/dist/lib/pkgconfig: pkg-config --cflags alps`
#INCL = -I${prefix}/include/coin-or
INCL += $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/home/james/dist/lib/pkgconfig: pkg-config --libs alps --static`
#LIBS = -L${exec_prefix}/lib -lAlps -lCoinUtils
LIBS += $(ADDLIBS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
