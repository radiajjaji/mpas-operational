#!/bin/bash
module load PrgEnv-intel
module swap intel intel-oneapi
module load cray-hdf5
module load cray-netcdf
module load cray-parallel-netcdf
module load craype-x86-rome
module load craype-hugepages2M
module load iobuf
module unload cuda

export NETCDF=$NETCDF_DIR
export WRFIO_NCD_LARGE_FILE_SUPPORT=1


cd /scratch/lus/arw/model/src/util/convert
./compile 
