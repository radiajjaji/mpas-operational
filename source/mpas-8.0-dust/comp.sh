#!/bin/bash
#SBATCH --error=/scratch/lus/arw/model/src/mpas/cce/8.0.t/compile.out
#SBATCH --output=/scratch/lus/arw/model/src/mpas/cce/8.0.t/compile.out
#SBATCH --job-name=CMP.MPS.CRAY
#SBATCH --partition=slm
#SBATCH --time=20:00:00
#SBATCH --ntasks=32
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL

set -euo pipefail
set -x

module purge

# --- Platform/target first ---
module load craype/2.7.32
module load craype-x86-rome
module load craype-network-ofi
module load libfabric/1.20.1

# --- Programming environment (this will pull the matching cce) ---
module load PrgEnv-cray/8.5.0

# --- MPI and math libs ---
module load cray-mpich/8.1.30
module load cray-libsci/24.07.0
module load cray-dsmml/0.3.0

# --- I/O libs (after PrgEnv + MPI are fixed) ---
module load cray-netcdf/4.9.0.13
module load cray-parallel-netcdf/1.12.3.13

# Avoid forcing MPAS to use an alternate I/O stack
unset PIO

module list

# Optional but recommended: makes runtime resolution consistent with your working ldd
export LD_LIBRARY_PATH=/opt/cray/pe/lib64:${LD_LIBRARY_PATH:-}

cd /scratch/lus/arw/model/src/mpas/cce/8.0.t/


# Clean + build init_atmosphere
make clean CORE=init_atmosphere
make -j 32 cray CORE=init_atmosphere CPP_EXTRA_FLAGS="-D DO_CHEMISTRY"

# Clean + build atmosphere
make clean CORE=atmosphere
make -j 1 cray CORE=atmosphere OPENMP=true CPP_EXTRA_FLAGS="-D DO_CHEMISTRY"
