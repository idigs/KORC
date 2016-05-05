#Makefile

UNAME	:=	$(shell uname)

SRC=./src
OBJ=./obj
MOD=./mod
BIN=./bin

ifeq ($(UNAME), Darwin)
	FC = mpif90

	FCDEFINEFLAGS =-DDOUBLE_PRECISION	-DHDF5_SINGLE_PRESICION
	FCFLAGS =-cpp   -J$(MOD) -fopenmp    -O3

	HDF5_INSTALL =../HDF5

	INCLUDES_HDF5	=-I$(HDF5_INSTALL)/include
	INCLUDES =
	EXTLIB	=-L$(HDF5_INSTALL)/lib
	LIBSHDF5	=$(EXTLIB)  $(HDF5_INSTALL)/lib/libhdf5_fortran.a   $(HDF5_INSTALL)/lib/libhdf5_hl.a    $(HDF5_INSTALL)/lib/libhdf5.a
	LIB_FLAGS =
	LIB	= -lm   -lz -ldl
	OMP_LIB_FLAGS =	-lgomp
endif

ifeq ($(UNAME), Linux)
	FC = mpif90

	FCDEFINEFLAGS =-DDOUBLE_PRECISION	-DHDF5_SINGLE_PRESICION
	FCFLAGS =-cpp	-J$(MOD)	-fopenmp	-O3

	HDF5_INSTALL =../HDF5

	INCLUDES_HDF5	=-I$(HDF5_INSTALL)/include
	INCLUDES =
	EXTLIB	=-L$(HDF5_INSTALL)/lib
	LIBSHDF5	=$(EXTLIB)	$(HDF5_INSTALL)/lib/libhdf5hl_fortran.a	$(HDF5_INSTALL)/lib/libhdf5_fortran.a	$(HDF5_INSTALL)/lib/libhdf5_hl.a	$(HDF5_INSTALL)/lib/libhdf5.a
	LIB_FLAGS =
	LIB	= -lm   -lz -ldl
	OMP_LIB_FLAGS =	-lgomp
endif

BINARIES = $(addprefix $(OBJ)/, korc_types.o    constants.o korc_hpc.o  external_subroutines.o	initialize.o    finalize.o	korc_units.o    emf.o	pic.o	korc_HDF5.o	main.o)

all: KORC-FO

$(OBJ)/korc_types.o: $(SRC)/korc_types.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/korc_types.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/constants.o: $(SRC)/constants.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/constants.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/external_subroutines.o: $(SRC)/external_subroutines.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/external_subroutines.f90  $(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/initialize.o: $(SRC)/initialize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/initialize.f90  $(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/finalize.o: $(SRC)/finalize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/finalize.f90  $(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_hpc.o: $(SRC)/korc_hpc.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/korc_hpc.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_HDF5.o: $(SRC)/korc_HDF5.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $(SRC)/korc_HDF5.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/pic.o: $(SRC)/pic.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/pic.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/emf.o: $(SRC)/emf.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/emf.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_units.o: $(SRC)/korc_units.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/korc_units.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/main.o: $(SRC)/main.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$(SRC)/main.f90		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@


KORC-FO: $(BINARIES)
	$(FC)	$(FCFLAGS)  $(BINARIES)	-o	$(BIN)/KORC-FO	$(LIB_FLAGS)	$(LIBSHDF5)	$(LIB)	$(OMP_LIB_FLAGS)


info:
	@echo	"COMPILING KORC-KO FOR: "$(UNAME)
	@echo	"STARTING COMPILATION..."

clean:
	@echo	"CLEANING..."
	rm -f $(BINARIES) $(MOD)/*.mod

