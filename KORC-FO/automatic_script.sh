#!/bin/bash
for it in 1 2 4 8 10
do
rm outputFiles/*h5
mpirun -np 1 -x OMP_NUM_THREADS=$it ./bin/KORC inputFiles/input_file.korc outputFiles/
sleep 1s
done
