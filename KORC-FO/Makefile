#Makefile

UNAME	:=	$(shell uname)

SRC=src
OBJ=obj
MOD=mod
BIN=bin

ifeq ($(UNAME), Darwin)
	FC = mpif90

	FCDEFINEFLAGS =-DDOUBLE_PRECISION	-DHDF5_DOUBLE_PRESICION
	FCFLAGS =-cpp   -J$(MOD) -fopenmp    -O2    -ffpe-trap=invalid

	HDF5_INSTALL =../HDF5
	PSPLINE_INSTALL=../pspline_lib

	INCLUDES_HDF5	=-I$(HDF5_INSTALL)/include
	INCLUDES =-I$(PSPLINE_INSTALL)/mod
	EXTLIB	=-L$(HDF5_INSTALL)/lib
	LIBSHDF5	=$(EXTLIB)  $(HDF5_INSTALL)/lib/libhdf5_fortran.a   $(HDF5_INSTALL)/lib/libhdf5_hl.a    $(HDF5_INSTALL)/lib/libhdf5.a
	LIB_FLAGS =-L$(PSPLINE_INSTALL)/lib
	LIB	= -lm   -lz -ldl	-lpspline
	OMP_LIB_FLAGS =	-lgomp
endif

ifeq ($(UNAME), Linux)
	#FC = mpiifor
	FC = ftn

	FCDEFINEFLAGS =-DDOUBLE_PRECISION	-DHDF5_DOUBLE_PRESICION	-UH5_DEBUG_API
	FCFLAGS =-cpp	-module $(MOD)	-openmp	-O2 #-mp

	HDF5_INSTALL =/opt/cray/hdf5/1.10.0/INTEL/15.0
	PSPLINE_INSTALL=../pspline_lib

	INCLUDES_HDF5	=-I$(HDF5_INSTALL)/include
	INCLUDES =-I$(PSPLINE_INSTALL)/mod #-I/opt/intel/impi/5.1.1.109/intel64/include
	EXTLIB	=-L$(HDF5_INSTALL)/lib
	LIBSHDF5	=$(EXTLIB)	$(HDF5_INSTALL)/lib/libhdf5hl_fortran.a	$(HDF5_INSTALL)/lib/libhdf5_fortran.a	$(HDF5_INSTALL)/lib/libhdf5_hl.a	$(HDF5_INSTALL)/lib/libhdf5.a
	LIB_FLAGS =-L$(PSPLINE_INSTALL)/lib #-L/opt/intel/impi/5.1.1.109/intel64/lib
	LIB	= -lm   -lz -ldl	-lpspline
	OMP_LIB_FLAGS =
endif

BINARIES = $(addprefix $(OBJ)/, korc_types.o    constants.o	korc_HDF5.o	korc_hpc.o  rnd_numbers.o	korc_fields.o	korc_coords.o	korc_interp.o	initialize.o    finalize.o	korc_units.o   korc_ppusher.o	main.o)


all: $(BIN)/KORC-FO

$(OBJ)/korc_types.o: $(SRC)/korc_types.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/constants.o: $(SRC)/constants.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/rnd_numbers.o: $(SRC)/rnd_numbers.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<  	$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/initialize.o: $(SRC)/initialize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<  	$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/finalize.o: $(SRC)/finalize.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<  	$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_hpc.o: $(SRC)/korc_hpc.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_HDF5.o: $(SRC)/korc_HDF5.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)    $<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_fields.o: $(SRC)/korc_fields.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_interp.o: $(SRC)/korc_interp.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_coords.o: $(SRC)/korc_coords.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_ppusher.o: $(SRC)/korc_ppusher.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/korc_units.o: $(SRC)/korc_units.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@

$(OBJ)/main.o: $(SRC)/main.f90
	$(FC)	-c	$(FCFLAGS)	$(FCDEFINEFLAGS)	$<		$(INCLUDES)	$(INCLUDES_HDF5)	$(LIB_FLAGS)	-o	$@


$(BIN)/KORC-FO: $(BINARIES)
	$(FC)	$(FCFLAGS)  $(BINARIES)	-o	$@	$(LIB_FLAGS)	$(LIBSHDF5)	$(LIB)	$(OMP_LIB_FLAGS)


info:
	@echo	"COMPILING KORC-KO FOR: "$(UNAME)
	@echo	"STARTING COMPILATION..."
	@echo $(MSG)

clean:
	@echo	"CLEANING..."
	rm -f $(MOD)/*.mod $(OBJ)/*.o

