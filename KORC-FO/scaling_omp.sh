#!/bin/bash
for it in `seq 10`
do
mpirun -np 1 -x OMP_NUM_THREADS=$it bin/KORC-FO inputFiles/input_file.korc outputFiles/scaling/
sleep 5s
done
