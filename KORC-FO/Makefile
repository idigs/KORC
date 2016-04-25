#Makefile

UNAME	:=	$(shell uname)

ifeq ($(UNAME), Darwin)
	FC = mpif90

	FCDEFINEFLAGS =
	FCFLAGS =   -J modules/   -ffixed-line-length-none    -fopenmp    -g  -fdefault-real-8    -fdefault-integer-8

	HDF_INSTALL =# ../HDF5

	INCLUDE_HDF5	=#	-I$(HDF_include)/INSTALL
	INCLUDES =
	EXTLIB	=# -L$(HDF_INSTALL)/lib
	LIBSHDF	=# $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 	
	LIB_FLAGS =
	LIB	= -lm
	OMP_LIB_FLAGS =	-lgomp
endif

ifeq ($(UNAME), Linux)
	FC = mpif90

	FCDEFINEFLAGS =#-DHDF5_FLOAT	-DHDF5_SINGLE_FILE
	FCFLAGS =-std=f95	-ffixed-line-length-none    -fopenmp    -g  -fdefault-real-8    -fdefault-integer-8

	HDF_INSTALL =# ../HDF5

	INCLUDE_HDF5	=#	-I$(HDF_include)/INSTALL
	INCLUDES =
	EXTLIB	=# -L$(HDF_INSTALL)/lib
	LIBSHDF	=# $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 	
	LIB_FLAGS =
	LIB	= -lm
	OMP_LIB_FLAGS =	-lgomp
endif

BINARIES = types.o	initialize.o	main_mpi.o	main.o

all: KORC-FO

main.o: sources/main.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	sources/main.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

main_mpi.o: sources/main_mpi.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	sources/main_mpi.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

types.o: sources/types.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	sources/types.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

#structures.o: sources/structures.f90
#	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	sources/structures.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

initialize.o: sources/initialize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	sources/initialize.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)

KORC-FO: $(BINARIES)
	$(FC)	$(FCFLAGS)  $(BINARIES)	-o	KORC-FO	$(LIB_FLAGS)	$(LIBSHDF)	$(LIB)	$(OMP_LIB_FLAGS)


info:
#	@echo	"Compiling the code for MacOS"
	@echo	$(UNAME)
	@echo	"COMPILING..."

clean:
	@echo	"CLEANING COMPILATION..."
	rm -f $(BINARIES) KORC-FO modules/*.mod

