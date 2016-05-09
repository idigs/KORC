#!/bin/bash
for it in 1 2 3 4 5 6 7 8 9 10 11 12
do
mpirun -np 1 -x OMP_NUM_THREADS=$it bin/KORC-FO inputFiles/input_file_scaling.korc outputFiles/
sleep 5s
done
