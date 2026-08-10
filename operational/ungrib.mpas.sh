#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/ungrib.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/ungrib.%a.out
#SBATCH --job-name=UNGRIB.MPAS --no-requeue
#SBATCH --ntasks=1
#SBATCH --partition=slm 
#SBATCH --export=ALL
#SBATCH --time=0:30:00
#SBATCH --exclusive
#SBATCH --distribution=block:block
#SBATCH --propagate=STACK,CORE
#SBATCH --array=0-1:3

# --- Setup Environment ---
>/scratch/lus/arw/model/log/ungrib.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

cd $TMPDIR

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
ulimit -a
ulimit -c 0


set -x

# --- Date Calculation ---
export CYCLE_RUN=$(expr $START_DATE | cut -c9-10)
export START_DATE=$($SMSDATE -12 $START_DATE)
export CYCLE=$(expr $START_DATE | cut -c9-10)
export YMD=$(expr $START_DATE | cut -c1-8)

[ ! -f $MODEL_RUN/r$CYCLE_RUN/OK.UNGRIB.MPAS.$START_DATE ] || \
/bin/rm $MODEL_RUN/r$CYCLE_RUN/OK.UNGRIB.MPAS.$START_DATE

# --- Filenames Setup ---
FHR=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
export LBC_FILE=gdas.t${CYCLE}z.pgrb2.0p25.f${FHR}
export CHEM_FILE=gefs.chem.t${CYCLE}z.a3d_0p50.f${FHR}.grib2
export MERGED_FILE=gdas_chem.t${CYCLE}z.f${FHR}.grib2

/bin/rm -v -f $LOCAL_LBC_DIR/lbc.$SLURM_ARRAY_TASK_ID

if [ $GET_LBC_LOCALLY = no ] ; then
    # 1. Download GFS (Atmosphere/GDAS)
    cat << EOF > FTP_GFS.SCR
    lcd $LOCAL_LBC_DIR
    open https://nomads.ncep.noaa.gov/
    cd /pub/data/nccf/com/gfs/prod/gdas.${YMD}/${CYCLE}/atmos
    repeat -d 10 --until-ok pget -c -n 8 gdas.t${CYCLE}z.pgrb2.0p25.f${FHR} -o ${LBC_FILE}
EOF
    lftp -f FTP_GFS.SCR

    # 2. Download GEFS-Chem (Aerosols)
    cat << EOF > FTP_CHEM.SCR
    lcd $LOCAL_LBC_DIR
    open https://nomads.ncep.noaa.gov/
    cd /pub/data/nccf/com/gens/prod/gefs.${YMD}/${CYCLE}/chem/pgrb2ap5
    repeat -d 10 --until-ok pget -c -n 8 ${CHEM_FILE}
EOF
    lftp -f FTP_CHEM.SCR

    # 3. Run the Python Merger
    echo "--> Running combine_gfs_chem.py for Hour ${FHR}..."
    $MODEL_ROOT/combine_gfs_chem.py \
        $LOCAL_LBC_DIR/$CHEM_FILE \
        $LOCAL_LBC_DIR/$LBC_FILE \
        $LOCAL_LBC_DIR/$MERGED_FILE

    if [ -f $LOCAL_LBC_DIR/$MERGED_FILE ]; then
        export LBC_TO_USE=$MERGED_FILE
    else
        echo "WARNING: Merger failed or output not found. Falling back to GFS-only."
        export LBC_TO_USE=$LBC_FILE
    fi
else
    # Logic for Local Copy
    if [ ! -f $GFS_LBC_DIR/r$CYCLE/$LBC_FILE ] ; then
       echo "Local GFS file missing"
       exit 1
    fi
    cp -v -p $GFS_LBC_DIR/r$CYCLE/$LBC_FILE $LOCAL_LBC_DIR/$LBC_FILE
    export LBC_TO_USE=$LBC_FILE
fi

# --- Configure Namelist ---
export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)

export EXE=$UNGRIB
export IOBUF_PARAMS='$MODEL_RUN/r$CYCLE/GFS*:verbose:count=12:size=32M,GRIBFILE.*:verbose:count=12:size=32M'

sed -e s"/YEAR/$ACTUAL_YEAR/"g \
    -e s"/MONTH/$ACTUAL_MONTH/"g \
    -e s"/DAY/$ACTUAL_DAY/"g \
    -e s"/HOUR/$ACTUAL_HOUR/"g \
    -e s"/DOMAIN/$WRF_DOMAIN/"g \
    -e s"/CYCLE/$CYCLE_RUN/"g \
    -e s"/VERSION/$WRF_VERSION/"g \
    -e s"/COMPILER/$COMPILER/"g \
    -e s"/LBC_TYPE/GFS/"g \
    -e s"/LBC_FREQUENCY/$((LBC_FREQUENCY*3600))/"g \
    $MODEL_ROOT/name/wrf/$WRF_VERSION/namelist.ungrib > namelist.wps

# --- Link Files and Run Ungrib ---
# We link the MERGED file (or fallback) to GRIBFILE.AAA
ln -sf $LOCAL_LBC_DIR/${LBC_TO_USE} GRIBFILE.AAA
ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFS Vtable
ln -sf $MODEL_LOG/ungrib.${FHR}.log ungrib.log

ldd $EXE
$STIME srun --cpu-bind=cores $EXE

# --- Check Completion ---
if [ $SLURM_ARRAY_TASK_ID -eq 0 ] ; then
   isuccess=0
   for range in $(seq $SLURM_ARRAY_TASK_MIN $LBC_FREQUENCY $SLURM_ARRAY_TASK_MAX ) ; do
       while true ; do
           grep "Successful completion of program ungrib.exe" $MODEL_LOG/ungrib.$(printf "%03d" $range).log
           if [ $? -eq 0 ] ; then
               isuccess=$((isuccess+1))
               touch $MODEL_RUN/r$CYCLE_RUN/OK.UNGRIB.MPAS.$START_DATE
               break
           else
               sleep 1
           fi
       done
   done
fi
