#!/bin/bash

#SBATCH --error=/home/oper/arw/model/src/post/intel/UPPV4.1/upp.err
#SBATCH --output=/home/oper/arw/model/src/post/intel/UPPV4.1/upp.out
#SBATCH --job-name=WRF.UPP --no-requeue
#SBATCH --nodes=1 --ntasks=8 --ntasks-per-node=8 --cpus-per-task=4
#SBATCH --partition=workq
#SBATCH --export=ALL
#SBATCH --time=00:10:00
#SBATCH --exclusive
#SBATCH --distribution=block:block
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL

#
>/home/oper/arw/model/src/post/intel/UPPV4.1/upp.err
>/home/oper/arw/model/src/post/intel/UPPV4.1/upp.out

set +x
export DATE=2021060603
export CYCLE=00
export COMPILER=intel
hostname
module restore PrgEnv-intel
module load cray-netcdf
set -x

cd $TMPDIR 

export domain=d01
export SMSDATE=/home/oper/dev/bin/smsdate
export MODEL_ROOT=/home/oper/arw/model
export WRF_DATA=/home/oper/arw/model/src/wrf/$COMPILER/4.2.2/WRF-4.2.2/run
export MODEL_RUN=/scratch/lus/tmp/arw/day/r$CYCLE

export tmmark=tm$CYCLE
export UNIPOST_HOME=$MODEL_ROOT/src/post/$COMPILER/UPPV4.1
export CRTMDIR=$UNIPOST_HOME/sorc/comlibs/crtm2/src/fix
export EXEC=$UNIPOST_HOME/exec/unipost.dev.exe

ln -sf $UNIPOST_HOME/parm/post_avblflds.xml .
ln -sf $UNIPOST_HOME/parm/post_avblflds_raphrrr.xml .
ln -sf $UNIPOST_HOME/parm/post_avblflds_comupp.xml .
ln -sf $UNIPOST_HOME/parm/params_grib2_tbl_new params_grib2_tbl_new 
ln -sf $UNIPOST_HOME/parm/nam_micro_lookup.dat .
ln -sf $UNIPOST_HOME/parm/hires_micro_lookup.dat .
ln -sf $UNIPOST_HOME/parm/postxconfig-NT-WRF.txt postxconfig-NT.txt
ln -sf $UNIPOST_HOME/parm/wrf_cntrl.parm wrf_cntrl.parm
ln -sf $UNIPOST_HOME/parm/wrf_cntrl.parm fort.14
ln -sf $UNIPOST_HOME/parm/optics_luts* .
ln -sf $CRTMDIR/EmisCoeff/IR_Water/Big_Endian/Nalli.IRwater.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM4.MWwater.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM5.MWwater.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM6.MWwater.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/IR_Land/SEcategory/Big_Endian/NPOESS.IRland.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/IR_Snow/SEcategory/Big_Endian/NPOESS.IRsnow.EmisCoeff.bin           ./
ln -sf $CRTMDIR/EmisCoeff/IR_Ice/SEcategory/Big_Endian/NPOESS.IRice.EmisCoeff.bin           ./
ln -sf $CRTMDIR/AerosolCoeff/Big_Endian/AerosolCoeff.bin     ./
ln -sf $CRTMDIR/CloudCoeff/Big_Endian/CloudCoeff.bin         ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_g11.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g11.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_g12.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g12.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_g13.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g13.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_g15.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g15.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_mt1r.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_mt1r.TauCoeff.bin
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_mt2.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_mt2.TauCoeff.bin
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/imgr_insat3d.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_insat3d.TauCoeff.bin
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/amsre_aqua.SpcCoeff.bin  ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/amsre_aqua.TauCoeff.bin  ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/tmi_trmm.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/tmi_trmm.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f13.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f13.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f14.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f14.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f15.SpcCoeff.bin    ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f15.TauCoeff.bin    ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f16.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f16.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f17.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f17.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f18.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f18.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f19.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f19.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f20.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f20.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/seviri_m10.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/seviri_m10.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/v.seviri_m10.SpcCoeff.bin   ./
ln -sf $CRTMDIR/TauCoeff/ODPS/Big_Endian/abi_gr.TauCoeff.bin   ./
ln -sf $CRTMDIR/SpcCoeff/Big_Endian/abi_gr.SpcCoeff.bin   ./
ln -sf griddef.out fort.110

aa=$(expr $DATE | cut -c1-4)
mm=$(expr $DATE | cut -c5-6)
jj=$(expr $DATE | cut -c7-8)
hh=$(expr $DATE | cut -c9-10)

cat > itag <<EOF
./wrfout_${domain}_${aa}-${mm}-${jj}_${hh}_00_00.nc
netcdf
${aa}-${mm}-${jj}_${hh}:00:00
NCAR
EOF

ln -sf $MODEL_RUN/wrfout_${domain}_${aa}-${mm}-${jj}_${hh}_00_00.nc .

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_ABORT_ON_ERROR=1
export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_HINTS="wrfout*:cb_nodes=8:striping_factor=8,WRF*:cb_nodes=8:striping_factor=8:cray_cb_write_lock_mode=2:cray_cb_nodes_multiplier=4"
export MPICH_RANK_REORDER_METHOD=0
export MPICH_PTL_MATCH_OFF=1
export MPICH_MPIIO_AGGREGATOR_PLACEMENT_DISPLAY=1

export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=134217728
export FOR_DISABLE_KMP_MALLOC=TRUE

export MPICH_RANK_REORDER_DISPLAY=1
export MPICH_OFI_VERBOSE=1
export MPICH_OFI_NIC_VERBOSE=1
export MPICH_OFI_NUM_NICS=1
export FI_OFI_RXM_TX_SIZE=128000
export FI_OFI_RXM_RX_SIZE=128000
export MPICH_COLL_OPT_OFF=MPI_Scatterv
export MPICH_COLL_SYNC=MPI_Gather


[ ! -f wrfout_${domain}_${aa}-${mm}-${jj}_${hh}_00_00.nc ] || \
time srun --cpu-bind=cores /home/oper/arw/model/src/post/$COMPILER/UPPV4.1/exec/unipost.dev.exe
mv WRFPRS* $WRF_DATA/wrfprs_${domain}_${aa}-${mm}-${jj}_${hh}_00_00.grb
