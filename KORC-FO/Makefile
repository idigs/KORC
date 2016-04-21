#Makefile

FC=mpif90

FLAGS=-std=c++0x
FC_DEFINE_FLAGS=-D
CCFLAGS =	$(CCXXFLAGS)	-fopenmp	-O3

HDF_INSTALL = ../HDF5

INCLUDE_HDF5	=	-I$(HDF_INSTALL)/include
INCLUDES =	-I../arma_libs/usr/include	-I../arma_libs/usr/lib	-I../arma_libs/usr/lib64
EXTLIB	= -L$(HDF_INSTALL)/lib
LIBSHDF	= $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 	
LIB_FLAGS =	-L../arma_libs/usr/lib -L../arma_libs/usr/lib64
LIB	=	-lz	-ldl	-lm	-lhdf5_cpp	-larmadillo
OMP_LIB_FLAGS =	-lgomp

BINARIES = main.o	main_mpi.o	types.o structures.o

all: KORC-FO

main.o: sources/main.cpp
	$(FC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/main.cpp		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

generalFunctions.o:	sources/generalFunctions.cpp	sources/generalFunctions.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/generalFunctions.cpp	$(INCLUDES)	$(LIB_FLAGS)

boundaryConditions.o:	sources/boundaryConditions.cpp	sources/boundaryConditions.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/boundaryConditions.cpp	$(INCLUDES)	$(LIB_FLAGS)

emf.o:	sources/emf.cpp	sources/emf.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/emf.cpp	$(INCLUDES)	$(LIB_FLAGS)

PIC.o:	sources/PIC.cpp	sources/PIC.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/PIC.cpp	$(INCLUDES)	$(LIB_FLAGS)

initialize.o:	sources/initialize.cpp	sources/initialize.h		sources/structures.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/initialize.cpp	$(INCLUDES)	$(LIB_FLAGS)

units.o:		sources/units.cpp	sources/units.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/units.cpp	$(INCLUDES)	$(LIB_FLAGS)

types.o:		sources/types.cpp	sources/types.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/types.cpp	$(INCLUDES)	$(LIB_FLAGS)

outputHDF5.o:		sources/outputHDF5.cpp	sources/outputHDF5.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/outputHDF5.cpp	$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	$(OMP_LIB_FLAGS)

quietStart.o:		sources/quietStart.cpp	sources/quietStart.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/quietStart.cpp	$(INCLUDES)	$(LIB_FLAGS)

randomStart.o:		sources/randomStart.cpp	sources/randomStart.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/randomStart.cpp	$(INCLUDES)	$(LIB_FLAGS)

alfvenic.o:		sources/alfvenic.cpp	sources/alfvenic.h
	$(CC)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	sources/alfvenic.cpp	$(INCLUDES)	$(LIB_FLAGS)

PROMETHEUS: $(BINARIES)
	$(CC)	$(CCFLAGS)  $(BINARIES)	-o	PROMETHEUS	$(LIB_FLAGS)	$(LIBSHDF)	$(LIB)	$(OMP_LIB_FLAGS)


info:
#	@echo	"Compiling the code for MacOS"
	@echo	$(UNAME)
	@echo	"COMPILING..."

clean:
	@echo	"CLEANING COMPILATION..."
	rm -f $(BINARIES)

