LIBRARY = libfruit.a
FCOMPILER = /opt/local/bin/gfortran-mp-8

all: ${OBJS}
	${FCOMPILER} -c *.f90
	ar cr ${LIBRARY} *.o

clean:
	find . -name "*.o" -or -name "*.mod" -or -name "${LIBRARY}" | xargs rm -rf

test:
	${FCOMPILER}


