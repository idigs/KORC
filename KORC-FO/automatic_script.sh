#!/bin/bash
for it in {60,50,40,30,20,10}
do
mpirun -np 1 -x OMP_NUM_THREADS=5 bin/KORC-FO inputFiles/SyntheticDiagnostic_pplane/input_file_${it}MeV.korc outputFiles/primary_generation/confinement/${it}MeV/
sleep 5s
done
