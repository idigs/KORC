#Makefile

UNAME:=$(shell uname)

SRC=src
OBJ=obj
MOD=mod
BIN=bin

ifeq ($(UNAME), Darwin)
        FC = mpif90

        FCDEFINEFLAGS =-DDOUBLE_PRECISION       -DHDF5_DOUBLE_PRESICION
        FCFLAGS =-cpp   -J$(MOD) -fopenmp    -O2

        PSPLINE_INSTALL=../pspline_lib

        HDF5_INSTALL =../HDF5

        INCLUDES_HDF5   =-I$(HDF5_INSTALL)/include
        INCLUDES =-I$(PSPLINE_INSTALL)/mod
        EXTLIB  =-L$(HDF5_INSTALL)/lib
        LIBSHDF5        =$(EXTLIB)  $(HDF5_INSTALL)/lib/libhdf5_fortran.a   $(HDF5_INSTALL)/lib/libhdf5_hl.a    $(HDF5_INSTALL)/lib/libhdf5.a
        LIB_FLAGS =-L$(PSPLINE_INSTALL)/lib
        LIB     = -lm   -lz -ldl        -lpspline
        OMP_LIB_FLAGS = -lgomp
endif #Darwin

ifeq ($(UNAME), Linux)

ifdef NERSC_HOST
        FC = ftn

ifdef PRG_ENV
ifeq ($(PRG_ENV),GNU)
        HDF5_INSTALL =../HDF5_GNU
        FCFLAGS =-cpp   -J$(MOD)        -fopenmp        -Ofast  #-O3    -ffast-math     -funroll-loops
        FCDEFINEFLAGS =-DDOUBLE_PRECISION       -DHDF5_DOUBLE_PRESICION
endif #GNU
ifeq ($(PRG_ENV),INTEL)
#        PSPLINE_INSTALL=../pspline_lib_INTEL
        HDF5_INSTALL =../HDF5_INTEL
        FCFLAGS =-cpp   -module $(MOD)  -qopenmp        -O2   #-g     -dynamic        -debug inline-debug-info
        FCDEFINEFLAGS =-DINTEL  -DDOUBLE_PRECISION      -DHDF5_DOUBLE_PRESICION
endif #INTEL
endif #PRG_ENV

ifeq ($(NERSC_HOST),cori)
        INCLUDES =$(PSPLINE)
        LIB_FLAGS =$(PSPLINE)
        LIB     = -lm   -lz -ldl        -lpspline
        OMP_LIB_FLAGS =
endif #cori

ifeq ($(NERSC_HOST),edison)
#       PSPLINE_INSTALL=../pspline_lib_GNU
#        INCLUDES =-I$(PSPLINE_INSTALL)/mod
#        LIB_FLAGS =-L$(PSPLINE_INSTALL)/lib
#        LIB     = -lm   -lz -ldl        -lpspline
        INCLUDES =$(PSPLINE)
        LIB_FLAGS =$(PSPLINE)
        LIB     = -lm   -lz -ldl -lpspline
        OMP_LIB_FLAGS =
endif #edison

else #NERSC_HOST
        FC = mpif90

        FCDEFINEFLAGS =-DDOUBLE_PRECISION       -DHDF5_DOUBLE_PRESICION -UH5_DEBUG_API
        FCFLAGS =-cpp   -J$(MOD)        -fopenmp        -O3

        HDF5_INSTALL =../HDF5
        PSPLINE_INSTALL=../pspline_lib

        INCLUDES =-I$(PSPLINE_INSTALL)/mod
        LIB_FLAGS =-L$(PSPLINE_INSTALL)/lib
endif #NERSC_HOST

endif # Linux

EXTLIB=-L$(HDF5_INSTALL)/lib
LIBSHDF5=$(EXTLIB)      $(HDF5_INSTALL)/lib/libhdf5hl_fortran.a $(HDF5_INSTALL)/lib/libhdf5_fortran.a   $(HDF5_INSTALL)/lib/libhdf5_hl.a        $(HDF5_INSTALL)/lib/libhdf5.a

INCLUDES_HDF5   =-I$(HDF5_INSTALL)/include
LIB     = -lm   -lz -ldl        -lpspline
OMP_LIB_FLAGS =

BINARIES = $(addprefix $(OBJ)/,korc_types.o     \
korc_constants.o        \
korc_HDF5.o     \
special_fun_modules.o   \
special_fun_subroutines.o       \
korc_synthetic_camera.o \
korc_collisions.o       \
korc_hpc.o      \
rnd_numbers.o   \
korc_coords.o   \
korc_interp.o   \
korc_fields.o   \
korc_avalanche.o        \
korc_experimental_pdf.o	\
hammersley_generator.o  \
korc_initialize.o       \
korc_finalize.o \
korc_units.o    \
korc_ppusher.o  \
main.o)



all: $(BIN)/KORC

$(OBJ)/%.o : $(SRC)/%.f90
	$(FC)   -c  $(FCFLAGS)  $(FCDEFINEFLAGS)    $<  $(INCLUDES) $(INCLUDES_HDF5)    $(LIB_FLAGS)    -o  $@

$(BIN)/KORC: $(BINARIES)
	$(FC)	$(FCFLAGS)  $(BINARIES)	-o	$@	$(LIB_FLAGS)	$(LIBSHDF5)	$(LIB)	$(OMP_LIB_FLAGS)


info:
	@echo	"COMPILING KORC-KO FOR: "$(UNAME)
	@echo	"COMPILING USING: "$(PRG_ENV)
ifdef NERSC_HOST
	@echo	"COMPILING AT: "$(NERSC_HOST)
else
	@echo	"COMPILING AT: "$(shell hostname)
endif
	@echo	"STARTING COMPILATION..."

clean:
	@echo	"CLEANING..."
	rm -f $(MOD)/*.mod	$(OBJ)/*.o	$(BIN)/KORC*

