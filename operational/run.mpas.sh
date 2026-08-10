#!/bin/bash

#SBATCH --job-name=MPAS.MODEL --requeue
#SBATCH --partition=opr 
#SBATCH --time=1:50:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL
#SBATCH --error=/scratch/lus/arw/model/log/mpas.out
#SBATCH --output=/scratch/lus/arw/model/log/mpas.out
#SBATCH --nodes=120 --ntasks=3720 --ntasks-per-node=31  --cpus-per-task=4

set -x

>/scratch/lus/arw/model/log/mpas.out
>/scratch/lus/arw/model/log/mpas.out
. /scratch/lus/arw/model/slurm/modelenv.sh

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
ulimit -a
ulimit -c 0
export KMP_AFFINITY=compact,verbose
export KMP_STACKSIZE=4G


cd $TMPDIR

export CYCLE_RUN=$(expr $START_DATE | cut -c9-10)
export START_DATE=$($SMSDATE -12 $START_DATE)
export START_YEAR=$(expr $START_DATE | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE | cut -c7-8)
export START_HOUR=$(expr $START_DATE | cut -c9-10)
export CYCLE=$(expr $START_DATE | cut -c9-10)
export DO_RESTART=false
export JEDI_RESTART_INTERVAL="${JEDI_RESTART_INTERVAL:-09:00:00}"

[ -f "$MODEL_RUN/r$CYCLE_RUN/OK.MPAS.$START_DATE" ] || [ "$(find "$MODEL_RUN/r$CYCLE_RUN" -maxdepth 1 -type f -name 'MPAS*.nc' | wc -l)" -eq 241 ] && exit 0

[ -f $MODEL_RUN/r$CYCLE_RUN/INIT.$START_YEAR-$START_MONTH-${START_DAY}_$START_HOUR.nc ] || exit 1

if [ $(ls -l $D/r$CYCLE_RUN/RSTR.*.nc | wc -l)  -gt 1 ] ; then
   export LAST_FILE=$(basename $(ls -rt1 $D/r$CYCLE_RUN/RSTR.*.nc | tail -1))
   export LAST_FILE_YEAR=$(expr  $LAST_FILE | cut -c6-9)
   export LAST_FILE_MONTH=$(expr $LAST_FILE | cut -c11-12)
   export LAST_FILE_DAY=$(expr   $LAST_FILE | cut -c14-15)
   export LAST_FILE_HOUR=$(expr  $LAST_FILE | cut -c17-18)

   date0=$(date -d "$START_YEAR/$START_MONTH/$START_DAY $START_HOUR:00:00" +%s)
   date1=$(date -d "$LAST_FILE_YEAR/$LAST_FILE_MONTH/$LAST_FILE_DAY $LAST_FILE_HOUR:00:00" +%s)

   last_range=$((date1-date0))
   last_range=$((last_range/10800))
   last_range=$((last_range*3))
   restart_range=$((MPAS_FORECAST_LENGTH_HOURS-last_range))

   export NDAYS=$((restart_range/24))
   export NHOURS=$((restart_range-24*NDAYS))
   export MPAS_FORECAST_LENGTH=${NDAYS}_${NHOURS}:00:00

   export START_YEAR=$LAST_FILE_YEAR
   export START_MONTH=$LAST_FILE_MONTH
   export START_DAY=$LAST_FILE_DAY
   export START_HOUR=$LAST_FILE_HOUR
   export DO_RESTART=true
fi


case $MPAS_RESOLUTION in
10)    XNN=x1; NUM=5898242;DT=60.0;DX=10000.0;UPWIND=0.0;PHOTDT=600 ;;
12)    XNN=x1; NUM=4096002;DT=72.0;DX=12000.0;UPWIND=0.0;PHOTDT=720 ;;
15)    XNN=x1; NUM=2621442;DT=90.0;DX=15000.0;UPWIND=0.0;PHOTDT=900 ;;
15_3)  XNN=x5; NUM=6488066;DT=20.0;DX=3000.0 ;UPWIND=0.0;PHOTDT=200 ;;
30_5)  XNN=x6; NUM=2819097;DT=36.0;DX=5000.0 ;UPWIND=0.0;PHOTDT=360 ;;
15_5)  XNN=x3; NUM=4763609;DT=36.0;DX=5000.0 ;UPWIND=0.0;PHOTDT=360 ;;
60_10) XNN=x6; NUM=999426; DT=60.0;DX=10000.0;UPWIND=0.5;PHOTDT=600 ;;
60_3)  XNN=x20;NUM=835586; DT=20.0;DX=3000.0 ;UPWIND=0.5;PHOTDT=200 ;;
46_12) XNN=x4; NUM=655362; DT=72.0;DX=12000.0;UPWIND=0.5;PHOTDT=720;;
60_15) XNN=x4; NUM=535554; DT=90.0;DX=15000.0;UPWIND=0.5;PHOTDT=900;;

*) exit 0;;
esac

export STRIDE=2
export STRIPE=8
export IOTASKS=$((SLURM_NTASKS/STRIPE))
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export IO_IN="pnetcdf,cdf5"
export IO_OUT="pnetcdf,cdf5"
export SSIZE=1048576
export MLTPSR=4
export MLTPSW=4
export MLTPW=6
export MLTPR=6
export MODE_LOCK=2
export MPICH_COLL_OPT_OFF=MPI_Scatterv
export MPICH_COLL_SYNC=MPI_Gather
export MPICH_RANK_REORDER_METHOD=0


export MPICH_MPIIO_HINTS="\
*/INIT.*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSR)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPR,\
*/MPAS.*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSW)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPW,\
*/RSTR.*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSW)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPW,\
*/DIAG.*:striping_factor=$STRIPE:striping_unit=$((SSIZE*MLTPSW)):cray_cb_write_lock_mode=$MODE_LOCK:cray_cb_nodes_multiplier=$MLTPW"

sed -e s"/START_DATE/$START_YEAR-$START_MONTH-${START_DAY}_${START_HOUR}:00:00/"g \
    -e s"/MPAS_FORECAST_LENGTH/$MPAS_FORECAST_LENGTH/"g \
    -e s"/UPWIND/$UPWIND/"g \
    -e s"/IOTASKS/$IOTASKS/"g \
    -e s"/STRIDE/$STRIDE/"g \
    -e s"/MPAS_RESOLUTION/$MPAS_RESOLUTION/"g \
    -e s"/PHOTDT/$PHOTDT/"g \
    -e s"/DT/$DT/"g \
    -e s"/DX/$DX/"g \
    -e s"/XNN/$XNN/"g \
    -e s"/NUM/$NUM/"g \
    -e s"/DO_RESTART/$DO_RESTART/"g \
    $MODEL_ROOT/name/mpas/$MPAS_VERSION/namelist.atmosphere > namelist.atmosphere

# Force +9h restart output for JEDI FGAT.
if grep -q "config_restart_interval" namelist.atmosphere; then
  sed -i "s/config_restart_interval *= *.*/config_restart_interval = '$JEDI_RESTART_INTERVAL'/" namelist.atmosphere
else
  echo "WARNING: config_restart_interval not found in namelist.atmosphere" >&2
fi

cat namelist.atmosphere

sed -e s"/CYCLE/$CYCLE_RUN/"g \
    -e s"/XNN/$XNN/"g \
    -e s"/NUM/$NUM/"g \
    -e s"/IO_IN/$IO_IN/"g \
    -e s"/IO_OUT/$IO_OUT/"g \
    $MODEL_ROOT/name/mpas/$MPAS_VERSION/streams.atmosphere > streams.atmosphere

cat streams.atmosphere

ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.diagnostics
ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.output
ln -sf $MODEL_ROOT/name/mpas/$MPAS_VERSION/stream_list.atmosphere.surface

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
ln -sf $MPAS_DIR/MP_THOMPSON_QRacrQS_DATA.DBL
ln -sf $MPAS_DIR/MP_THOMPSON_QRacrQG_DATA.DBL
ln -sf $MPAS_DIR/MP_THOMPSON_QIautQS_DATA.DBL
ln -sf $MPAS_DIR/MP_THOMPSON_freezeH2O_DATA.DBL

ln -sf $MPAS_DIR/build_tables
ln -sf $MPAS_DIR/default_inputs

/bin/rm -f $MODEL_LOG/log.atmosphere.0000.out

ln -sf $MPAS_DIR/atmosphere_model mpas.exe
ln -sf $MODEL_LOG/log.atmosphere.0000.out log.atmosphere.0000.out


$STIME srun -v -l --cpu-bind=cores --cpus-per-task=$OMP_NUM_THREADS mpas.exe

grep "Finished running the atmosphere core" $MODEL_LOG/log.atmosphere.0000.out



if [ $? -eq 0 ] ; then

   # After successful MPAS completion, move the +9h restart for JEDI FGAT.
   # During the run, RSTR files remain in r$CYCLE_RUN for crash recovery.
   ASSIM_DATE=$($SMSDATE 9 $START_DATE)
   ASSIM_YEAR=$(expr  $ASSIM_DATE | cut -c1-4)
   ASSIM_MONTH=$(expr $ASSIM_DATE | cut -c5-6)
   ASSIM_DAY=$(expr   $ASSIM_DATE | cut -c7-8)
   ASSIM_HOUR=$(expr  $ASSIM_DATE | cut -c9-10)

   ASSIM_DIR="$D/assim"
   mkdir -p "$ASSIM_DIR"

   SRC_RESTART="$D/r$CYCLE_RUN/RSTR.${ASSIM_YEAR}-${ASSIM_MONTH}-${ASSIM_DAY}_${ASSIM_HOUR}.nc"
   DST_RESTART="$ASSIM_DIR/restart.${ASSIM_YEAR}-${ASSIM_MONTH}-${ASSIM_DAY}_${ASSIM_HOUR}:00:00.nc"

   if [ -s "$SRC_RESTART" ]; then
      rm -f "$DST_RESTART"
      mv -f "$SRC_RESTART" "$DST_RESTART"
      echo "JEDI FGAT restart moved: $SRC_RESTART -> $DST_RESTART"
   else
      echo "WARNING: expected +9h restart not found: $SRC_RESTART" >&2
   fi

   touch $MODEL_RUN/r$CYCLE_RUN/OK.MPAS.$START_DATE
   exit 0
else
   cp $MODEL_LOG/log.atmosphere.0000.out $HOME/model/fail/.
   cp $MODEL_LOG/mpas.out $HOME/model/fail/.
   exit 1
fi

