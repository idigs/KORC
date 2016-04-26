#Makefile

UNAME	:=	$(shell uname)

SRC=./src
OBJ=./obj
MOD=./mod
BIN=./bin

ifeq ($(UNAME), Darwin)
	FC = mpif90

	FCDEFINEFLAGS = -DWITH_MPI
	FCFLAGS =   -J$(MOD) -fopenmp    -g  -cpp

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
	FCFLAGS =-J$(MOD)	-fopenmp    -g -fopenmp    -g  -cpp

	HDF_INSTALL =# ../HDF5

	INCLUDE_HDF5	=#	-I$(HDF_include)/INSTALL
	INCLUDES =
	EXTLIB	=# -L$(HDF_INSTALL)/lib
	LIBSHDF	=# $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 	
	LIB_FLAGS =
	LIB	= -lm
	OMP_LIB_FLAGS =	-lgomp
endif

BINARIES = $(addprefix $(OBJ)/, types.o  initialize.o    finalize.o  main_mpi.o  main.o)

all: KORC-FO

$(OBJ)/types.o: $(SRC)/types.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/types.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/initialize.o: $(SRC)/initialize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/initialize.f90  $(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/finalize.o: $(SRC)/finalize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/finalize.f90  $(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/main_mpi.o: $(SRC)/main_mpi.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/main_mpi.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/main.o: $(SRC)/main.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/main.f90		$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS)	-o	$@


KORC-FO: $(BINARIES)
	$(FC)	$(FCFLAGS)  $(BINARIES)	-o	$(BIN)/KORC-FO	$(LIB_FLAGS)	$(LIBSHDF)	$(LIB)	$(OMP_LIB_FLAGS)


info:
#	@echo	"Compiling the code for MacOS"
	@echo	$(UNAME)
	@echo	"COMPILING..."

clean:
	@echo	"CLEANING COMPILATION..."
	rm -f $(BINARIES) $(MOD)/*.mod

