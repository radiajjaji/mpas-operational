#!/bin/bash

#SBATCH --job-name=MPAS.INIT --no-requeue
#SBATCH --partition=opr
#SBATCH --export=ALL
#SBATCH --time=0:10:00
#SBATCH --exclusive
#SBATCH --error=/scratch/lus/arw/model/log/init.mpas.out
#SBATCH --output=/scratch/lus/arw/model/log/init.mpas.out
#SBATCH --nodes=16 --ntasks=256 --ntasks-per-node=16

>/scratch/lus/arw/model/log/init.mpas.out
>/scratch/lus/arw/model/log/init.mpas.out
. /scratch/lus/arw/model/slurm/modelenv.sh

cd $TMPDIR

set +x

module purge
unset LD_LIBRARY_PATH
unset LIBRARY_PATH
unset CPATH
unset CMAKE_PREFIX_PATH
unset PKG_CONFIG_PATH

if [ $COMPILER_MPAS = intel ] ; then
   module load PrgEnv-intel/8.5.0
fi

if [ $COMPILER_MPAS = cce ] ; then
   module load PrgEnv-cray/8.5.0
fi

if [ $COMPILER_MPAS = aocc ] ; then
   module load PrgEnv-aocc/8.5.0
fi

if [ $COMPILER_MPAS = nvhpc ] ; then
   module load PrgEnv-nvhpc/8.5.0
fi

module load cray-netcdf/4.9.0.13
module load cray-parallel-netcdf/1.12.3.13
module load craype-x86-rome

if [ $USE_UCX = yes ] ; then
   module swap craype-network-ofi craype-network-ucx
   module swap cray-mpich cray-mpich-ucx
fi

module list

ulimit -s unlimited
set -x

export CYCLE_RUN=$(expr $START_DATE | cut -c9-10)
export START_DATE=$($SMSDATE -12 $START_DATE)
export START_YEAR=$(expr $START_DATE | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE | cut -c7-8)
export START_HOUR=$(expr $START_DATE | cut -c9-10)
export CYCLE=$(expr $START_DATE | cut -c9-10)
[ ! -f $MODEL_RUN/r$CYCLE_RUN/OK.INIT.MPAS.$START_DATE ] || \
/bin/rm $MODEL_RUN/r$CYCLE_RUN/OK.INIT.MPAS.$START_DATE

STRIPE=8
SSIZE=1048576

export MPICH_MPIIO_HINTS="*/GFS*:cb_nodes=32,*/INIT.*:cb_nodes=$STRIPE:striping_factor=$STRIPE:striping_unit=$SSIZE"
export MPICH_COLL_OPT_OFF=MPI_Scatterv
export MPICH_COLL_SYNC=MPI_Gather
export MPICH_RANK_REORDER_METHOD=0

case $MPAS_RESOLUTION in
10)    XNN=x1; NUM=5898242; LEVS=56 ;;
12)    XNN=x1; NUM=4096002; LEVS=56 ;;
15)    XNN=x1; NUM=2621442; LEVS=56 ;;
15_3)  XNN=x5; NUM=8060930; LEVS=56 ;;
15_5)  XNN=x3; NUM=4763609; LEVS=56 ;;
30_5)  XNN=x6; NUM=2819097; LEVS=56 ;;
60_10) XNN=x6; NUM=999426;  LEVS=56 ;;
60_3)  XNN=x20;NUM=835586;  LEVS=56 ;;
46_12) XNN=x4; NUM=655362;  LEVS=56;;
60_15) XNN=x4; NUM=535554;  LEVS=56;;
*) exit 0;;
esac

export IO_IN="pnetcdf,cdf5"
export IO_OUT="pnetcdf,cdf5"

export STRIDE=10
export IOTASKS=16
export NPROC=$SLURM_NTASKS
export PPN=$SLURM_NTASKS_PER_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export KMP_STACKSIZE=2G
export KMP_AFFINITY=disabled

sed -e s"/START_DATE/$START_YEAR-$START_MONTH-${START_DAY}_${START_HOUR}:00:00/"g \
    -e s"/CYCLE/$CYCLE_RUN/"g \
    -e s"/MPAS_LEVS/$MPAS_LEVS/"g \
    -e s"/MPAS_RESOLUTION/$MPAS_RESOLUTION/"g \
    -e s"/GFS_LEVS/$GFS_LEVS/"g \
    -e s"/IOTASKS/$IOTASKS/"g \
    -e s"/STRIDE/$STRIDE/"g \
    -e s"/XNN/$XNN/"g \
    -e s"/NUM/$NUM/"g \
    $MODEL_ROOT/name/mpas/$MPAS_VERSION/namelist.init_atmosphere > namelist.init_atmosphere

cat namelist.init_atmosphere

sed -e s"/IO_IN/$IO_IN/"g \
    -e s"/MPAS_RESOLUTION/$MPAS_RESOLUTION/"g \
    -e s"/XNN/$XNN/"g \
    -e s"/NUM/$NUM/"g \
    -e s"/IO_OUT/$IO_OUT/"g \
    -e s"/CYCLE/$CYCLE_RUN/"g \
    $MODEL_ROOT/name/mpas/$MPAS_VERSION/streams.init_atmosphere > streams.init_atmosphere
    
cat streams.init_atmosphere

ln -sf $MPAS_DIR/init_atmosphere_model init.exe
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/CAM_ABS_DATA.DBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/CAM_AEROPT_DATA.DBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/GENPARM.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/LANDUSE.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/OZONE_DAT.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/OZONE_LAT.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/OZONE_PLEV.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/RRTMG_LW_DATA
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/RRTMG_LW_DATA.DBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/RRTMG_SW_DATA
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/RRTMG_SW_DATA.DBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/SOILPARM.TBL
ln -sf $MPAS_DIR/src/core_atmosphere/physics/physics_wrf/files/VEGPARM.TBL
ln -sf $MPAS_DIR/build_tables
ln -sf $MPAS_DIR/default_inputs

ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.diagnostics
ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.output
ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.surface

[ ! -f $MODEL_RUN/r$CYCLE_RUN/INIT.$START_YEAR-$START_MONTH-${START_DAY}_${START_HOUR}.nc ] || \
/bin/rm -f $MODEL_RUN/r$CYCLE_RUN/INIT.$START_YEAR-$START_MONTH-${START_DAY}_${START_HOUR}.nc
/bin/rm -f $MODEL_LOG/log.init_atmosphere.0000.out

export SLURM_CPU_BIND=verbose
ln -sf $MODEL_LOG/log.init_atmosphere.0000.out log.init_atmosphere.0000.out

ldd init.exe

$STIME srun -l \
   --export=ALL \
   --propagate=STACK,CORE \
   --cpu_bind=core \
   init.exe

# Make another try ... slurm + cray shasta are sometimes crazy

if [ $? -ne 0 ] ; then
   $STIME srun -l \
   --export=ALL \
   --propagate=STACK,CORE \
   --cpu_bind=core \
   --cpus-per-task=$OMP_NUM_THREADS \
   init.exe
fi

grep "Finished running the init_atmosphere core" $MODEL_LOG/log.init_atmosphere.0000.out
if [ $? -eq 0 ] ; then
   touch $MODEL_RUN/r$CYCLE_RUN/OK.INIT.MPAS.$START_DATE
fi
