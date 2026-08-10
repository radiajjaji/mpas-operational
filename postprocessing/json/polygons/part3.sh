#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --job-name=JSN.MPS.3
#SBATCH --ntasks=8
#SBATCH --partition=opr --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

# Set environment variables
set -x
. /scratch/lus/arw/model/slurm/modelenv.sh

export ATTRNAME='DN'
export WIND_THIN_FACTOR=2

cd $TMPDIR

echo "Running on :" $(hostname)
export CYCLE=$(expr $START_DATE | cut -c9-10)

export START_YEAR=$(expr $START_DATE | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE | cut -c7-8)
export START_HOUR=$(expr $START_DATE | cut -c9-10)
export RUN_DATE=$(expr $START_DATE | cut -c1-8)

export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)
export GRBFILE=$D/r$CYCLE/mpsprs.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.grb1

[ -f $GRBFILE ] || exit 0

$GRIB2CTL -verf $GRBFILE > mpsprs.ctl
$GRIBMAP -big -i mpsprs.ctl

# Define all parameters

export parameters=("DBZ" "BRTMP" "APCP6H" "GP700C" "GP700" "DIR700" "GP500" "GP500C" "DIR500" "SPDJET" "DIRJET" "CAPE" "CIN" "LINDEX" "VERVEL")

if [ $SLURM_ARRAY_TASK_ID -ge 0 -a $((SLURM_ARRAY_TASK_ID % 6)) -eq 0 ] ; then
   INC=6
   if [ $SLURM_ARRAY_TASK_ID -eq 0 ] ; then
	   INC=0
   fi
   export DATE_6H=$($SMSDATE -$INC $ACTUAL_DATE)
   export YEAR_6H=$(expr $DATE_6H | cut -c1-4)
   export MONTH_6H=$(expr $DATE_6H | cut -c5-6)
   export DAY_6H=$(expr $DATE_6H | cut -c7-8)
   export HOUR_6H=$(expr $DATE_6H | cut -c9-10)
   export GRBFILE_6H=$D/r$CYCLE/mpsprs.$YEAR_6H-$MONTH_6H-${DAY_6H}_${HOUR_6H}_00_00.grb1
   $GRIB2CTL -verf $GRBFILE > mpsprs.ctl
   $GRIBMAP -big -i mpsprs.ctl
   $GRIB2CTL -verf $GRBFILE_6H > mpsprsp.ctl
   $GRIBMAP -big -i mpsprsp.ctl
else
   export parameters=("DBZ" "BRTMP" "GP700C" "GP700" "DIR700" "GP500" "GP500C" "DIR500" "SPDJET" "DIRJET" "CAPE" "CIN" "LINDEX" "VERVEL")
fi

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

# Process each parameter

for param in "${parameters[@]}"; do

  if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
     [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

  case $param in
    "APCP6H")
      export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
      time $GDALCALC -A $GRBFILE --A_band=$index -B $GRBFILE_6H --B_band=$index --outfile=$param.tiff --calc="A-B"
      time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME
        ;;
    "DBZ")
      export index=$($WGRIB $GRBFILE | grep "REFC" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "BRTMP")
      export index=$($WGRIB $GRBFILE | grep "BRTMP" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "GP700"|"GP500")
      export index=$($WGRIB $GRBFILE | grep "HGT" | grep "${param:2:3} mb" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "SPDJET")
      export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "250 mb" | awk -F ":" '{print $1}')
      export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "250 mb" | awk -F ":" '{print $1}')
      time $GDALCALC -A $GRBFILE --A_band $index_u -B $GRBFILE --B_band $index_v --outfile $param.tiff \
           --calc="(1.994*sqrt(A*A+B*B))*((1.994*sqrt(A*A+B)) > 50)"
      time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "CAPE"|"CIN")
      export index=$($WGRIB $GRBFILE | grep "$param" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "LINDEX")
      export index=$($WGRIB $GRBFILE | grep "4LFTX" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;
    "VERVEL")
      export index=$($WGRIB $GRBFILE | grep "VVEL" | awk -F ":" '{print $1}')
      time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
      ;;

    "DIR700"|"DIR500")
      export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "${param:3:3} mb" | awk -F ":" '{print $1}')
      export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "${param:3:3} mb" | awk -F ":" '{print $1}')
      $WGRIB -d $index_u -grib -o u.grb $GRBFILE
      $WGRIB -d $index_v -grib -o v.grb $GRBFILE
      $CDO -nint -mulc,1.994 -sqrt -add -sqr u.grb -sqr v.grb spd.grb
      $CDO -nint -addc,180. -mulc,57.3 -atan2 -mulc,-1 u.grb -mulc,-1 v.grb dir.grb
      $GRIB2CTL -verf spd.grb > spd.ctl
      $GRIBMAP -big -i spd.ctl
      $GRIB2CTL -verf dir.grb > dir.ctl
      $GRIBMAP -big -i dir.ctl
      cat << EOF > $param.gs
'open spd.ctl'
'open dir.ctl'
'set gxout shp'
'set shp -pt $param'
'd skip(1000*ugrd${param:3:3}mb.1+ugrd${param:3:3}mb.2,$WIND_THIN_FACTOR)'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp
      sed -i s"/GRID_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.geojson
      ;;

    "DIRJET")
      export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "250 mb" | awk -F ":" '{print $1}')
      export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "250 mb" | awk -F ":" '{print $1}')
      $WGRIB -d $index_u -grib -o u.grb $GRBFILE
      $WGRIB -d $index_v -grib -o v.grb $GRBFILE
      $CDO -nint -mulc,1.994 -sqrt -add -sqr u.grb -sqr v.grb spd.grb
      $CDO -nint -addc,180. -mulc,57.3 -atan2 -mulc,-1 u.grb -mulc,-1 v.grb dir.grb
      $GRIB2CTL -verf spd.grb > spd.ctl
      $GRIB2CTL -verf dir.grb > dir.ctl
      $GRIBMAP -big -i spd.ctl
      $GRIBMAP -big -i dir.ctl
      cat << EOF > $param.gs
'open spd.ctl'
'open dir.ctl'
'set gxout shp'
'set shp -pt $param'
'd skip(1000*ugrd250mb.1+ugrd250mb.2,$WIND_THIN_FACTOR)'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      sed -i s"/GRID_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.geojson
      ;;
    "GP700C"|"GP500C")
      cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set lev ${param:2:3}'
'set cint 20'
'd hgtprs'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Line.geojson $param.shp ; sed -i s"/CNTR_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Line.geojson
      /bin/rm $param.shp $param.gs

      cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -pt $param -fmt 10 1'
'set lev ${param:2:3}'
'd skip(hgtprs,20)'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Point.geojson $param.shp ; sed -i s"/GRID_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Point.geojson
      /bin/rm $param.shp $param.gs

      if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
         time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Line.geojson
         time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.Point.mbtiles $OUTPUTPATH/$param.Point.geojson
         time $TILEJOIN -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Point.mbtiles
	 
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
      ;;
  esac
  fi
  if [ "$param" != "GP700C" ] && [ "$param" != "GP500C" ] && [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
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
