#!/bin/tcsh
#PBS -N ad:dsmr_SimTests
#PBS -l nodes=6:c18c:ppn=12
#PBS -l walltime=24:00:00
##PBS -q alpha
#PBS -j oe

cd $PBS_O_WORKDIR

mpirun --map-by node -np 72 python-mpi /sciclone/home00/geogdan/geoMatch_testing/SimTests/simInit.py
