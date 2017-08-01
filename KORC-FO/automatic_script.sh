#!/bin/bash
for it in $(seq 1 1 10)
do
mpirun -np 1 -x OMP_NUM_THREADS=10 ./bin/KORC inputFiles/input_file.korc /media/l8c/FantomHD/DIII-D/KORC_data/test/
sleep 1s
done
