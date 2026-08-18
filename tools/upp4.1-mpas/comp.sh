#!/bin/bash

module load PrgEnv-intel
module switch intel intel-oneapi
module load cray-netcdf
module load cray-mpich

export NETCDF=$NETCDF_DIR
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

cd /scratch/lus/arw/model/src/post/intel/UPPV4.1
./compile

