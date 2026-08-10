# Initialize Module environment

#source /usr/local/Modules/init/bash
source /scratch/lus/$USER/.bash_profile > /dev/null 2>&1

# Start Date

export START_DATE=2026081012

# Model Environment variables

export MPAS_WITH_DUST=yes
export USE_WRF_RESTART_FACLITY=yes
export COMPILER=intel
export COMPILER_MPAS=cce
export COMPILER_WRF=cce
export COMPILER_WRFDA=intel
export COMPILER_JEDI=intel
export WRF_VERSION=4.7.1
export MPAS_RESOLUTION=12
export MPAS_LEVS=56
export WRFVAR_CRTM_VERSION=2.3.0
export GSI_CRTM_VERSION=2.3.0
export WRFDA_VERSION=4.7.1
export GSI_VERSION=3.7

export LBC_TYPE=GFS
export LBC_TYPE=MPAS
export SIZE_METGRID_MPAS=7083982280
export SIZE_METGRID_GFS=7083983428
export SIZE_UNGRIB_MPAS=5754291504
export SIZE_UNGRIB_GFS=5754291504
export SIZE_UNGRIB_VAR="SIZE_UNGRIB_${LBC_TYPE}"
export SIZE_UNGRIB=${!SIZE_UNGRIB_VAR}
export SIZE_METGRID_VAR="SIZE_METGRID_${LBC_TYPE}"
export SIZE_METGRID=${!SIZE_METGRID_VAR}

export DIRECTLY_FROM_MPAS=no
export NUM_METGRID_LEVELS=42

if [ $DIRECTLY_FROM_MPAS = yes ] ; then
   export NUM_METGRID_LEVELS=57
fi

if [ $MPAS_WITH_DUST = yes ] ; then
   export MPAS_VERSION=8.0.v
   export DUST=".dust"
   export MPASOUT_FILE_SIZE=24000000000
else
   export MPAS_VERSION=8.2.2
   export DUST=""
   export MPASOUT_FILE_SIZE=11000000000
fi

export NDAYS=7
export FORECAST_LENGTH=168
export LBC_FREQUENCY=3
export GFS_MAX_FORECAST=180
export ARCHIVE_FREQUENCY=3
export GET_LBC_LOCALLY=no
export FORECAST_LENGTH_ASSIM=9
export LBC_FREQUENCY_ASSIM=3
export BUFR=1
export CRTM=1
export EM_CORE=1
export WRF_EM_CORE=1
export NMM_CORE=0
export COAMPS_CORE=0
export WRF_NMM_CORE=0
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export WRF_CHEM=1
export WRF_KPP=1
export WRF_DOMAIN=AFAD
export WLON=-180.
export ELON=180.
export SLAT=-90.
export NLAT=90.
export MPAS_DOMAIN=GLOBAL
export MPAS_FORECAST_LENGTH=10_00:00:00
export MPAS_FORECAST_LENGTH_HOURS=240
export GFS_LEVS=42

export STAT_NMC=yes
export CHEM_IN_OPT=0
export HYBRID_OPT=2
export USE_THETA_M=1
export IO_API=PNETCDF
export CV_OPTIONS=3
export FGAT=yes
export BLEND=yes
export FIRST_GUESS=MPAS
export WIND_THIN=5
export OGR_GEOMETRY_ACCEPT_UNCLOSED_RING=YES
export TILESERVER=maptiler/tileserver-gl
export TILESERVER=tileserver-gl-light

# Model Different Paths

export SOFT_ROOT=/scratch/lus/dev
export MODEL_ROOT=/scratch/lus/arw/model
export ASCII_ROOT=/scratch/lus/arw/model/src/ascii
export OBS_SRC=/scratch/lus/arw/model/src/obs
export MODEL_RUN=/scratch/lus/tmp/arw/day/
export WORKDIR=/scratch/lus/tmp
export ARCHIVE=/scratch/lus/archive/
export INSTALLROOT=$MODEL_ROOT/wps
export EXT_DATAROOT=$MODEL_ROOT/wps/domains/AFAD
export MOAD_DATAROOT=$MODEL_ROOT/wps/domains/AFAD
export JASPERLIB=$SOFT_ROOT/lib
export JASPERINC=$SOFT_ROOT/include
export RTTOV=$MODEL_ROOT/var/rttov
export MODEL_LOG=$MODEL_ROOT/log
export LOG=$MODEL_ROOT/log
export LASTLOG=$MODEL_ROOT/lastlog
export WRF_DIR=$MODEL_ROOT/src/wrf/$COMPILER_WRF/$WRF_VERSION/WRF
export VAR_DIR=$MODEL_ROOT/src/wrf/$COMPILER_WRFDA/$WRFDA_VERSION/WRFDA/var
export WPS_DIR=$MODEL_ROOT/src/wrf/$COMPILER_WRF/$WRF_VERSION/WPS
export GSI_DIR=$MODEL_ROOT/src/gsi/$COMPILER/$GSI_VERSION
export GSI_FIX_DIR=$MODEL_ROOT/src/gsi/$COMPILER/$GSI_VERSION/fix
export CRTM_DIR=$MODEL_ROOT/src/crtm
export WPS_GEOG=$MODEL_ROOT/wps/WPS_GEOG
export WRF_DOM=$MODEL_ROOT/wps/domains/$WRF_DOMAIN/geogrid/modis/$WRF_VERSION
export MPAS_DIR=$MODEL_ROOT/src/mpas/$COMPILER_MPAS/$MPAS_VERSION
export MESH_DIR=$MODEL_ROOT/src/mpas/meshes/$MPAS_RESOLUTION/
export UNIPOST_DIR=$MODEL_ROOT/src/post/$COMPILER/UPPV4.1
export UNIPOST_CRTM_DIR=$UNIPOST_DIR/sorc/comlibs/crtm2/src/fix
export POST_OUT_FORMAT=grib2
export GSI_CRTM=$MODEL_ROOT/src/crtm/$GSI_CRTM_VERSION
export GFS_LBC_DIR=/scratch/lus/tmp/arw/day/gfs
export MPAS_LBC_DIR=/scratch/lus/tmp/arw/day
export LOCAL_LBC_DIR=/scratch/lus/tmp/arw/day/lbc
export ARCHIVE_WRF=/scratch/lus/archive/wrf
export ARCHIVE_MPAS=/scratch/lus/archive/mps
export ARCHIVE_GFS=/scratch/lus/archive/gfs
export REMOTE_ARCHIVE=/path/to/remote/archive
export GRAPH_WRF=/scratch/lus/data/WRF
export REMOTE_GRAPH_MPAS=/opt/MPS
export REMOTE_GRAPH_WRF=/opt/WRF
export REMOTE_GRAPH_GFS=/path/to/remote/GFS
export REMOTE_GRAPH_AIGFS=/path/to/remote/AIGFS
export REMOTE_GRAPH_AIFS=/path/to/remote/AIFS
export REMOTE_GRAPH_PANGU=/path/to/remote/PANGU
export REMOTE_GRAPH_GCAST=/path/to/remote/GCAST
export GRAPH_MPAS=/scratch/lus/data/MPS
export GRAPH_GFS=/scratch/lus/data/GFS
export GRAPH_AIGFS=/scratch/lus/data/AIGFS
export GRAPH_AIFS=/scratch/lus/data/AIFS
export GRAPH_PANGU=/scratch/lus/data/PANGU
export GRAPH_GCAST=/scratch/lus/data/GCAST
export GRAPH_ECM=/scratch/lus/data/ECM
export GRAPH_ICN=/scratch/lus/data/ICN
export GRIB_LOC=/scratch/lus/data/LOC
export SMS_ROOT=/scratch/lus/tmp/arw/day/sms/
export T=$WORKDIR
export R=$MODEL_ROOT/log
export B=$WORKDIR/arw/day/lbc
export D=$WORKDIR/arw/day
export O=$ARCHIVE
export DAD=$MODEL_RUN/assim
export ECMWF_HOST="https://data.ecmwf.int/forecasts"
export ECMWF_DIR="ifs/0p25/oper"

# External Resources Sites

export NCEP=ftpprd.ncep.noaa.gov
export NCEP_DIR=/pub/data/nccf/com/obsproc/prod/
export NCEP_SITE1=ftp://ftpprd.ncep.noaa.gov/
export NCEP_SITE2=https://nomads.ncep.noaa.gov
export NCEP_SITE=$NCEP_SITE2

# Model utils

export JAVA_ROOT=/usr/lib64/jvm/java
export JAVA_HOME=/usr/lib64/jvm/java
export JAVA_BINDIR=/usr/lib64/jvm/java/bin
export PB2LTR=$MODEL_ROOT/src/obs/converter/pb2ltr.exe
export LTR2NDG=$MODEL_ROOT/src/obs/converter/RT_fdda_reformat_obsnud.pl
export SMSDATE=/scratch/lus/dev/bin/smsdate
export JSMIN=/scratch/lus/dev/bin/jsmin
export JSMIN="jq -c . "
export SQUEUE=/usr/bin/squeue
export SCANCEL=/usr/bin/scancel
export SBATCH=/usr/bin/sbatch
export STIME="/usr/bin/time -v"
export LFTP=/scratch/lus/dev/bin/lftp
export NCL=/scratch/lus/dev/bin/ncl
export NCATTED=/scratch/lus/dev/bin/ncatted
export COPYGB="/scratch/lus/dev/bin/copygb.exe"
export PCOPYGB="/scratch/lus/dev/bin/pcopygb.py"
export CDO=/scratch/lus/dev/bin/cdo
export WGRIB2=/scratch/lus/dev/bin/wgrib2
export WGRIB=/scratch/lus/dev/bin/wgrib
export GRB1TO2="/scratch/lus/dev/bin/grb1to2.pl "
export GRIB2CTL=/scratch/lus/dev/bin/grib2ctl.pl
export GBLAV2CTL=/scratch/lus/arw/bin/gblav2ctl.pl
export G2CTL=/scratch/lus/dev/bin/g2ctl.pl
export ALTG2CTL=/scratch/lus/dev/bin/alt_g2ctl
export GRIBMAP=/scratch/lus/dev/opengrads/gribmap
export ALTGRIBMAP=/scratch/lus/dev/bin/alt_gmp
export PB2NC=$MODEL_ROOT/src/met.new/bin/pb2nc
export POINTSTAT=$MODEL_ROOT/src/met.new/bin/point_stat
export GRIB2JSON=/scratch/lus/dev/bin/grib2json
export CNVGRIB=/scratch/lus/dev/bin/cnvgrib
export PCNVGRIB=/scratch/lus/dev/bin/pcnvgrib.py
export GDALPOLYGON=/scratch/lus/dev/bin/gdal_polygonize.py
export OGR2OGR=/scratch/lus/dev/bin/ogr2ogr
export OGRMERGE=/scratch/lus/dev/bin/ogrmerge.py
export LABELS=/scratch/lus/dev/bin/lab.py
export GDALWARP="/scratch/lus/dev/bin/gdalwarp -t_srs '+proj=latlong +datum=WGS84' "
export GDALWARP="/scratch/lus/dev/bin/gdalwarp -t_srs EPSG:4326 "
export MAX_ZOOM_TILES=7
export OBS2IODA=/scratch/lus/arw/model/src/obs2ioda/obs2ioda.x
export TILEJOIN="/scratch/lus/dev/bin/tile-join -pk "
export GDALCONT=/scratch/lus/dev/bin/gdal_contour
export GDAL_TRANSLATE=/scratch/lus/dev/bin/gdal_translate
export GDALBAND="/scratch/lus/dev/bin/gdal_contour -p"
export ISOBANDS="/scratch/lus/dev/bin/isobands.matplot_lib"
export RASTER_TO_GEOJSON="/scratch/lus/dev/bin/raster_to_geojson "
export GDALCALC=/scratch/lus/dev/bin/gdal_calc.py
export OGR2OGR=/scratch/lus/dev/bin/ogr2ogr
export OGRMERGE=/scratch/lus/dev/bin/ogrmerge.py
export DOCKER=/usr/bin/docker
export GDALPERF=" -multi --config GDAL_CACHEMAX 512 --config GDAL_NUM_THREADS ALL_CPUS "
export GRADS="/scratch/lus/dev/opengrads/grads"
export MBUTIL="/scratch/lus/dev/bin/mb-util --image_format=pbf "
export MONTHS=("jan" "feb" "mar" "apr" "may" "jun" "jul" "aug" "sep" "oct" "nov" "dec")

export UNGRIB=$MODEL_ROOT/src/wrf/intel/4.2.2/WPS/orig/WPS-master/ungrib/src/ungrib.exe
export METGRID_MPAS=$MODEL_ROOT/src/util/metgrid_mpas/metgrid/src/metgrid.exe

export TIPPECANOE_MAX_THREADS=4
export TPCANOE="/scratch/lus/dev/bin/tippecanoe \
  --no-feature-limit \
  --temporary-directory=/dev/shm \
  --extend-zooms-if-still-dropping \
  --no-tiny-polygon-reduction \
  --no-simplification-of-shared-nodes \
  -P -f -z$MAX_ZOOM_TILES"

# CRAY MPICH Optimization

export MPICH_ENV_DISPLAY=0
export MPICH_VERSION_DISPLAY=0
export PMI_VERSION_DISPLAY=0
export MPICH_ABORT_ON_ERROR=1
export MPICH_RANK_REORDER_DISPLAY=0
export MPICH_MPIIO_HINTS_DISPLAY=0
export MPICH_MPIIO_AGGREGATOR_PLACEMENT_DISPLAY=0

# Memory Cray Interconnect Optimization

export USE_UCX=yes
export USE_OLD_MPICH=no
export MPICH_OFI_USE_PROVIDER="verbs;ofi_rxm"
export FI_OFI_RXM_TX_SIZE=128000
export FI_OFI_RXM_RX_SIZE=128000

# Stack management

ulimit -s unlimited
ulimit -a
ulimit -c 0

# Machines

export WEB_HOST=WEB_HOSTNAME_OR_IP
export WEB=WEB_HOSTNAME_OR_IP
export CLM1=CLUSTER_HOST_1
export CLM2=CLUSTER_HOST_2
export HP1=POSTPROCESS_HOST_1
export HP2=POSTPROCESS_HOST_2
export HP3=POSTPROCESS_HOST_3
export WEB_DIR=/path/to/web/root
export NWPOUT_HOST=NWP_OUTPUT_HOST
export NWPOUT_USER=nwpout
