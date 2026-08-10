#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/metgrid.mpas.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/metgrid.mpas.%a.out
#SBATCH --job-name=MTGRD.MPAS
#SBATCH --nodes=1 --ntasks=32
#SBATCH --partition=opr  --no-requeue
#SBATCH --export=ALL
#SBATCH --time=0:15:00
#SBATCH --array=12-180:3

>>/scratch/lus/arw/model/log/metgrid.mpas.$SLURM_ARRAY_TASK_ID.out
>>/scratch/lus/arw/model/log/metgrid.mpas.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

cd $TMPDIR

set +x

module purge
unset LD_LIBRARY_PATH
unset LIBRARY_PATH
unset CPATH
unset CMAKE_PREFIX_PATH
unset PKG_CONFIG_PATH
module load craype/2.7.32

if [ $COMPILER = intel ] ; then
   module load PrgEnv-intel/8.5.0
fi

if [ $COMPILER = cce ] ; then
   module load PrgEnv-cray/8.5.0
fi

if [ $COMPILER = aocc ] ; then
   module load PrgEnv-aocc/8.5.0
fi

if [ $COMPILER = nvhpc ] ; then
   module load PrgEnv-nvhpc/8.5.0
fi

module load cray-mpich/8.1.30
module load cray-libsci/24.07.0
module load cray-dsmml/0.3.0
module load cray-hdf5/1.14.3.1
module load cray-netcdf/4.9.0.13
module load cray-parallel-netcdf/1.12.3.13
module load libfabric/1.20.1
module load craype-x86-rome
module load craype-network-ofi
module load atp

module list


ulimit -s unlimited
ulimit -c 0
ulimit -a unlimited
set -x

echo "Running on :" $(hostname)

export CYCLE=$(expr $START_DATE | cut -c9-10)
[ ! -f $MODEL_RUN/r$CYCLE/OK.METGRID.MPAS.$START_DATE ] || \
/bin/rm $MODEL_RUN/r$CYCLE/OK.METGRID.MPAS.$START_DATE

export START_DATE=$($SMSDATE -12 $START_DATE)
export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)

if [ $USE_UCX = yes ] ; then
   module swap craype-network-ofi craype-network-ucx
   module swap cray-mpich cray-mpich-ucx
fi

export EXE=$METGRID_MPAS


export SSIZE=1048576
export MLTPSR=4
export MLTPSW=4
export MLTPSW=4
export MLTPW=6
export MLTPR=6
export MODE_LOCK=2
export MPICH_COLL_OPT_OFF=MPI_Scatterv
export MPICH_COLL_SYNC=MPI_Gather
export MPICH_RANK_REORDER_METHOD=0

export MPICH_MPIIO_HINTS="\
*/MPAS*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSR)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPR,\
*/met_em*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSW)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPW"

export MPICH_PTL_MATCH_OFF=1

export MPICH_COLL_OPT_OFF=MPI_Scatterv
export MPICH_COLL_SYNC=MPI_Gather

export ATP_ENABLED=0
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G

case $MPAS_RESOLUTION in
10)    XNN=x1; NUM=5898242;DT=60.0;DX=10000.0;UPWIND=0.0;DUST=$DUST ;;
12)    XNN=x1; NUM=4096002;DT=72.0;DX=12000.0;UPWIND=0.0;DUST=$DUST ;;
15)    XNN=x1; NUM=2621442;DT=90.0;DX=15000.0;UPWIND=0.0;DUST=$DUST ;;
15_3)  XNN=x5; NUM=8060930;DT=20.0;DX=3000.0 ;UPWIND=0.0;DUST=$DUST ;;
15_5)  XNN=x3; NUM=4763609;DT=36.0;DX=5000.0 ;UPWIND=0.0;DUST=$DUST ;;
30_5)  XNN=x6; NUM=2819097;DT=36.0;DX=5000.0 ;UPWIND=0.0;DUST=$DUST ;;
60_10) XNN=x6; NUM=999426; DT=60.0;DX=10000.0;UPWIND=0.5;DUST=$DUST ;;
60_3)  XNN=x20;NUM=835586; DT=20.0;DX=3000.0 ;UPWIND=0.5;DUST=$DUST ;;
46_12) XNN=x4; NUM=655362; DT=72.0;DX=12000.0;UPWIND=0.5;DUST=$DUST;;
60_15) XNN=x4; NUM=535554; DT=90.0;DX=15000.0;UPWIND=0.5;DUST=$DUST;;

*) exit 0;;
esac

export MPAS_GRID=$XNN.$NUM.static$DUST.nc

sed -e s"/YEAR/$ACTUAL_YEAR/"g \
    -e s"/MONTH/$ACTUAL_MONTH/"g \
    -e s"/DAY/$ACTUAL_DAY/"g \
    -e s"/HOUR/$ACTUAL_HOUR/"g \
    -e s"/DOMAIN/$WRF_DOMAIN/"g \
    -e s"/MPAS_GRID/$MPAS_GRID/"g \
    -e s"/CYCLE/$CYCLE/"g \
    -e s"/VERSION/$WRF_VERSION/"g \
    -e s"/COMPILER/$COMPILER/"g \
    -e s"/LBC_FREQUENCY/$((LBC_FREQUENCY*3600))/"g \
    $MODEL_ROOT/name/wrf/$WRF_VERSION/namelist.metgrid.mpas > namelist.wps

[ ! -f $MODEL_RUN/r$CYCLE/met_em.d01.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.nc ] || exit 0

ln -sf $MODEL_LOG/metgrid.mpas.$(printf "%03d" $SLURM_ARRAY_TASK_ID).log metgrid.log.0000
ln -sf $MESH_DIR/$MPAS_GRID .

if [ -f $MODEL_RUN/r$CYCLE/MPAS.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_$ACTUAL_HOUR.nc \
     -a  1$(stat -c %s $MODEL_RUN/r$CYCLE/MPAS.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_$ACTUAL_HOUR.nc) -ge 1$MPASOUT_FILE_SIZE ] ; then
   ln -sf $MODEL_RUN/r$CYCLE/MPAS.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_$ACTUAL_HOUR.nc .
else
   sleep 5
fi

cat namelist.wps

$STIME srun -v -n $SLURM_NTASKS -N $SLURM_NNODES  --cpu-bind=cores $EXE
