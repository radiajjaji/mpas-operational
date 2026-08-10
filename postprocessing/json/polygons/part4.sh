#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --job-name=JSN.MPS.4
#SBATCH --ntasks=8
#SBATCH --partition=opr  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part4.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

export ATTRNAME='DN'

cd $TMPDIR
export CYCLE=$(expr $START_DATE | cut -c9-10)
set +x

module load PrgEnv-intel
module load cray-netcdf
module load cray-hdf5
module load craype-x86-rome
module load cray-libsci

ulimit -s unlimited
ulimit -c 0
ulimit -a unlimited
set -x

echo "Running on :" $(hostname)

export START_YEAR=$(expr $START_DATE  | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE   | cut -c7-8)
export START_HOUR=$(expr $START_DATE  | cut -c9-10)
export RUN_DATE=$(expr $START_DATE    | cut -c1-8)

export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)
export GRBFILE=$D/r$CYCLE/mpsprs.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.grb1

[ -f $GRBFILE ] || exit 0

$GRIB2CTL -verf $GRBFILE > mpsprs.ctl
$GRIBMAP -big -i mpsprs.ctl


# Input And Output Directories
export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

export params=("APCP24H" "RH700" "OLR" "PWATER" "VOR500" "CEIL" "DUST" "DVIS")

if [ $SLURM_ARRAY_TASK_ID -ge 0 -a $(( ($SLURM_ARRAY_TASK_ID + $CYCLE) % 24 )) -eq 0 ] ; then
   INC=24
   if [ $SLURM_ARRAY_TASK_ID -eq 0 ] ; then
           INC=$CYCLE
   fi
   export DATE_24H=$($SMSDATE -$INC $ACTUAL_DATE)
   export YEAR_24H=$(expr $DATE_24H | cut -c1-4)
   export MONTH_24H=$(expr $DATE_24H | cut -c5-6)
   export DAY_24H=$(expr $DATE_24H | cut -c7-8)
   export HOUR_24H=$(expr $DATE_24H | cut -c9-10)
   export GRBFILE_24H=$D/r$CYCLE/mpsprs.$YEAR_24H-$MONTH_24H-${DAY_24H}_${HOUR_24H}_00_00.grb1
   $GRIB2CTL -verf $GRBFILE_24H > mpsprsp.ctl
   $GRIBMAP -big -i mpsprsp.ctl
else
   export params=("RH700" "OLR" "PWATER" "VOR500" "CEIL" "DUST" "DVIS")
fi

# Process each parameter

for param in "${params[@]}"; do
  if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
     [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

  case $param in
       "APCP24H")
          export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
          time $GDALCALC -A $GRBFILE --A_band=$index -B $GRBFILE_24H --B_band=$index --outfile=$param.tiff --calc="A-B"
	  time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME
          ;;
       "RH700")
          export index=$($WGRIB $GRBFILE | grep "RH" | grep "${param:2:3} mb" | awk -F ":" '{print $1}')
	  time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
          ;;
       "OLR")
          export index=$($WGRIB $GRBFILE | grep "ULWRF" | grep "top" | awk -F ":" '{print $1}')
	  time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
          ;;
        "PWATER")
            export index=$($WGRIB $GRBFILE | grep "PWAT" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "CEIL")
            export index=$($WGRIB $GRBFILE | grep "HGT" | grep "cld base" | awk -F ":" '{print $1}')
            time $GDALCALC -A $GRBFILE --A_band=$index --calc="(A/1000.)"  --outfile=$param.tiff
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "DUST")
            cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout geotiff'
'set geotiff $param '
'define thresh = 0.15 - 0.02 * sltypsfc + 0.0035 * soilm0_200cm + 0.15 * vegsfc/100. + 0.05 * log(1 + sfcrsfc)'
'define dust = 100 * sltypsfc * pow(fricvsfc,3) * pow((1 - (thresh / fricvsfc)),2) * exp(-vegsfc/100.) * exp(-soilm0_200cm) * exp(-sfcrsfc)'
'define dust = maskout(dust,dust)'
'd dust'
'quit'
EOF
            time $GRADS -lbc $param.gs
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "DVIS")
            cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout geotiff'
'set geotiff $param '
'define thresh = 0.15 - 0.02 * sltypsfc + 0.0035 * soilm0_200cm + 0.15 * vegsfc/100. + 0.05 * log(1 + sfcrsfc)'
'define dust = 100 * sltypsfc * pow(fricvsfc,3) * pow((1 - (thresh / fricvsfc)),2) * exp(-vegsfc/100.) * exp(-soilm0_200cm) * exp(-sfcrsfc)'
'define dust = maskout(dust,dust)'
'define visib = 10 * exp(-0.1*dust)'
'd visib'
'quit'
EOF
            time $GRADS -lbc $param.gs
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "VOR500")
            cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout geotiff'
'set geotiff $param '
'set lev 500'
'd 10000*hcurl(ugrdprs,vgrdprs)'
'quit'
EOF
            time $GRADS -lbc $param.gs
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
    esac
  fi

  if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
     time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
     
CNT=$(sqlite3 "$OUTPUTPATH/$param.mbtiles" "select count(*) from tiles;" 2>/dev/null || echo 0)

if [ "$CNT" -le 0 ]; then
  echo "SKIP $param: empty tiles table => $OUTPUTPATH/$param.mbtiles"
else
  DEST="$REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/"
  rsync -avz --delay-updates --partial-dir=".rsync-partial" \
    "$OUTPUTPATH/$param.mbtiles" \
    "$USER@$WEB:$DEST"
fi
  fi
done
