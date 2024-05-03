#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_PSPLINE=ON \
      -DUSE_FIO=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-DHDF5_DOUBLE_PRESICION" \
    ../

make -j VERBOSE=1
