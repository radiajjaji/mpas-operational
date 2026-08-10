#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/upp.mpas.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/upp.mpas.%a.out
#SBATCH --job-name=UPP.MPAS
#SBATCH --nodes=1 --ntasks-per-node=8 --cpus-per-task=14
#SBATCH --partition=opr --exclusive
#SBATCH --time=0:20:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

>>/scratch/lus/arw/model/log/upp.mpas.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh


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

export OMP_PROC_BIND=close
export OMP_PLACES=threads
export KMP_STACKSIZE=4G
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

if [ $USE_UCX = yes ] ; then
   module swap craype-network-ofi craype-network-ucx
   module swap cray-mpich cray-mpich-ucx
fi

set -x

echo "Running on :" $(hostname)

cd $TMPDIR 

export CYCLE=$(expr $START_DATE | cut -c9-10)
export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)

export tmmark=tm$CYCLE
export EXEC=$UNIPOST_DIR/exec/unipost.mpas.exe
: "${PCNVGRIB:?PCNVGRIB is not set}"

ln -sf $MODEL_ROOT/etc/input_stations_global input_stations
ln -sf $UNIPOST_DIR/parm/post_avblflds.xml .
ln -sf $UNIPOST_DIR/parm/post_avblflds_raphrrr.xml .
ln -sf $UNIPOST_DIR/parm/post_avblflds_comupp.xml .
ln -sf $UNIPOST_DIR/parm/params_grib2_tbl_new params_grib2_tbl_new 
ln -sf $UNIPOST_DIR/parm/nam_micro_lookup.dat .
ln -sf $UNIPOST_DIR/parm/hires_micro_lookup.dat .
ln -sf $UNIPOST_DIR/parm/postxconfig.mpas.txt postxconfig-NT.txt
ln -sf $UNIPOST_DIR/parm/optics_luts* .
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/IR_Water/Big_Endian/Nalli.IRwater.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/MW_Water/Big_Endian/FASTEM4.MWwater.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/MW_Water/Big_Endian/FASTEM5.MWwater.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/MW_Water/Big_Endian/FASTEM6.MWwater.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/IR_Land/SEcategory/Big_Endian/NPOESS.IRland.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/IR_Snow/SEcategory/Big_Endian/NPOESS.IRsnow.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/EmisCoeff/IR_Ice/SEcategory/Big_Endian/NPOESS.IRice.EmisCoeff.bin           ./
ln -sf $UNIPOST_CRTM_DIR/AerosolCoeff/Big_Endian/AerosolCoeff.bin     ./
ln -sf $UNIPOST_CRTM_DIR/CloudCoeff/Big_Endian/CloudCoeff.bin         ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_g11.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_g11.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_g12.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_g12.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_g13.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_g13.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_g15.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_g15.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_mt1r.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_mt1r.TauCoeff.bin
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_mt2.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_mt2.TauCoeff.bin
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/imgr_insat3d.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/imgr_insat3d.TauCoeff.bin
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/amsre_aqua.SpcCoeff.bin  ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/amsre_aqua.TauCoeff.bin  ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/tmi_trmm.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/tmi_trmm.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmi_f13.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmi_f13.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmi_f14.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmi_f14.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmi_f15.SpcCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmi_f15.TauCoeff.bin    ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmis_f16.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmis_f16.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmis_f17.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmis_f17.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmis_f18.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmis_f18.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmis_f19.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmis_f19.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/ssmis_f20.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/ssmis_f20.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/seviri_m10.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/seviri_m10.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/v.seviri_m10.SpcCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/TauCoeff/ODPS/Big_Endian/abi_gr.TauCoeff.bin   ./
ln -sf $UNIPOST_CRTM_DIR/SpcCoeff/Big_Endian/abi_gr.SpcCoeff.bin   ./
ln -sf griddef.out fort.110
cp  $UNIPOST_DIR/parm/wrf_cntrl.mpas.parm wrf_cntrl.parm
cp  $UNIPOST_DIR/parm/wrf_cntrl.mpas.parm fort.14

aa=$(expr $ACTUAL_DATE | cut -c1-4)
mm=$(expr $ACTUAL_DATE | cut -c5-6)
jj=$(expr $ACTUAL_DATE | cut -c7-8)
hh=$(expr $ACTUAL_DATE | cut -c9-10)

export JULYR=$aa
export JULDAY=$(date -d "$aa-$mm-${jj}" +%j)

export START_YEAR=$(expr $START_DATE | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE | cut -c7-8)
export START_HOUR=$(expr $START_DATE | cut -c9-10)

if [ $SLURM_ARRAY_TASK_ID -eq 0 ] ; then
	aap=$aa
	mmp=$mm
        jjp=$jj
        hhp=$hh
else
        aap=$($SMSDATE -1 $ACTUAL_DATE | cut -c1-4)
        mmp=$($SMSDATE -1 $ACTUAL_DATE | cut -c5-6)
        jjp=$($SMSDATE -1 $ACTUAL_DATE | cut -c7-8)
        hhp=$($SMSDATE -1 $ACTUAL_DATE | cut -c9-10)
fi

cat > itag_grib1 <<EOF
./mpsout.d01.${aa}-${mm}-${jj}_${hh}_00_00.nc
netcdf
${aa}-${mm}-${jj}_${hh}_00_00
NCAR
EOF

cat > itag_grib2 <<EOF
./mpsout.d01.${aa}-${mm}-${jj}_${hh}_00_00.nc
netcdf
grib2
${aa}-${mm}-${jj}_${hh}_00_00
NCAR
EOF

export MPICH_MPIIO_HINTS="*/mpsout.d01.*:cb_nodes=8:striping_factor=6,*/MPASPR*:cb_nodes=8:striping_factor=8,*/mpsprs.*:cb_nodes=8:striping_factor=6"
export MPICH_RANK_REORDER_METHOD=0
export MPICH_PTL_MATCH_OFF=1

TARGET="$MODEL_RUN/r$CYCLE/mpsout.d01.${aa}-${mm}-${jj}_${hh}_00_00.nc"

for ((i=0; i<180; i++)); do
    [[ -f $TARGET && $(du -sk "$TARGET" | awk '{print $1}') -gt 20000000 ]] && {
        ln -sf "$TARGET" .
        echo "File ready: $TARGET"
        break
    }
    sleep 30
done

((i == 180)) && {
    echo "File not found or too small after 1 hour. Exiting."
    exit 0
}

# Generate Grib1 files
 
if [ ! -s $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1 ] ; then
   [ -d $GRIB_LOC/$START_DATE ] || mkdir $GRIB_LOC/$START_DATE
   cp itag_grib1 itag
   $STIME srun -v -N $SLURM_NNODES -n $SLURM_NTASKS --cpus-per-task=$OMP_NUM_THREADS --cpu-bind=cores $EXEC
   mv MPASPR* $GRIB_LOC/$START_DATE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1
   ln -sf $GRIB_LOC/$START_DATE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1 $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1
   $GRIB2CTL -verf $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1 > $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.ctl
   $GRIBMAP -big -i $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.ctl
fi
 
# Generate Grib 2 file

if [ ! -s $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb ] ; then
   $STIME $PCNVGRIB -g12 -nv -p40 -j ${SLURM_NTASKS:-1} --tmpdir "$TMPDIR" \
          $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb1 \
          $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
fi


# Generate  json wind data for javascript animations
 
mkdir -p $GRAPH_MPAS/windjson/$START_DATE/SRFC/

if [ ! -s $GRAPH_MPAS/windjson/$START_DATE/SRFC/${aa}${mm}${jj}${hh}.json ] ; then
   $STIME $WGRIB2 -match "10 m above ground" -GRIB wind.grb $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
   $STIME $CDO remapbil,$MODEL_ROOT/etc/gfs.griddes wind.grb wind.reduced.grb
   $STIME $GRIB2JSON --names --data --fp wind --fs 103 --fv 10.0 wind.reduced.grb > wind.json
   $JSMIN wind.json > $GRAPH_MPAS/windjson/$START_DATE/SRFC/${aa}${mm}${jj}${hh}.json
   /bin/rm wind.grb wind.reduced.grb wind.json
fi

mkdir -p $GRAPH_MPAS/windjson/$START_DATE/850/

if [ ! -s $GRAPH_MPAS/windjson/$START_DATE/850/${aa}${mm}${jj}${hh}.json ] ; then
   $STIME $WGRIB2 -match "UGRD:850|VGRD:850" -GRIB wind.grb $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
   $STIME $CDO remapbil,$MODEL_ROOT/etc/gfs.griddes wind.grb wind.reduced.grb
   $STIME $GRIB2JSON --data  wind.reduced.grb > wind.json
   $JSMIN wind.json > $GRAPH_MPAS/windjson/$START_DATE/850/${aa}${mm}${jj}${hh}.json
   /bin/rm wind.grb wind.reduced.grb wind.json
fi

mkdir -p $GRAPH_MPAS/windjson/$START_DATE/500/

if [ ! -s $GRAPH_MPAS/windjson/$START_DATE/500/${aa}${mm}${jj}${hh}.json ] ; then
   $STIME $WGRIB2 -match "UGRD:500|VGRD:500" -GRIB wind.grb $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
   $STIME $CDO remapbil,$MODEL_ROOT/etc/gfs.griddes wind.grb wind.reduced.grb
   $STIME $GRIB2JSON --data wind.reduced.grb > wind.json
   $JSMIN wind.json > $GRAPH_MPAS/windjson/$START_DATE/500/${aa}${mm}${jj}${hh}.json
   /bin/rm wind.grb wind.reduced.grb wind.json
fi

mkdir -p $GRAPH_MPAS/windjson/$START_DATE/700/

if [ ! -s $GRAPH_MPAS/windjson/$START_DATE/700/${aa}${mm}${jj}${hh}.json ] ; then
   $STIME $WGRIB2 -match "UGRD:700|VGRD:700" -GRIB wind.grb $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
   $STIME $CDO remapbil,$MODEL_ROOT/etc/gfs.griddes wind.grb wind.reduced.grb
   $STIME $GRIB2JSON --data wind.reduced.grb > wind.json
   $JSMIN wind.json > $GRAPH_MPAS/windjson/$START_DATE/700/${aa}${mm}${jj}${hh}.json
   /bin/rm wind.grb wind.reduced.grb wind.json
fi

mkdir -p $GRAPH_MPAS/windjson/$START_DATE/250/

if [ ! -s $GRAPH_MPAS/windjson/$START_DATE/250/${aa}${mm}${jj}${hh}.json ] ; then
   $STIME $WGRIB2 -match "UGRD:250|VGRD:250" -GRIB wind.grb $MODEL_RUN/r$CYCLE/mpsprs.${aa}-${mm}-${jj}_${hh}_00_00.grb
   $STIME $CDO remapbil,$MODEL_ROOT/etc/gfs.griddes wind.grb wind.reduced.grb
   $STIME $GRIB2JSON --data  wind.reduced.grb > wind.json
   $JSMIN wind.json > $GRAPH_MPAS/windjson/$START_DATE/250/${aa}${mm}${jj}${hh}.json
   /bin/rm wind.grb wind.reduced.grb wind.json
fi
