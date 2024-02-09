#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_

export CC=/usr/bin/gcc-12
export CXX=/usr/bin/g++-12
export FC=/usr/bin/gfortran-12

cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_OMP=OFF \
      -DUSE_PSPLINE=OFF \
      -DUSE_FIO=OFF \
      -DCORI_DIR=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-O3 -DHDF5_DOUBLE_PRESICION -fopenmp" \
      -DCMAKE_C_FLAGS="-O3 -fopenmp"  \
      -DCMAKE_CXX_FLAGS="-O3 -fopenmp" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-g -ffpe-trap=invalid,zero,overflow -fbacktrace -Werror" \
      -DCMAKE_C_FLAGS_DEBUG="-g" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g"

make VERBOSE=1
