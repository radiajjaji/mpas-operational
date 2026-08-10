#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --job-name=JSN.MPS.2
#SBATCH --ntasks=8
#SBATCH --partition=opr  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part2.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

export WIND_THIN_FACTOR=2
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

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

# Parameters to process
declare -a params=("GUST" "LCLD" "MCLD" "HCLD" "CBCLD" "VISIB" "LWRAD" "TROP" "FRZH" "TMP850" "TMP500" "GP850C" "GP850" "DIR850")

for param in "${params[@]}" ; do
  if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
     [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson


    case $param in
        "GUST")
            export index=$($WGRIB $GRBFILE | grep "GUST" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
            ;;
        "LCLD")
            export index=$($WGRIB $GRBFILE | grep "LCDC" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "MCLD")
            export index=$($WGRIB $GRBFILE | grep "MCDC" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "HCLD")
            export index=$($WGRIB $GRBFILE | grep "HCDC" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "CBCLD")
            export index_b=$($WGRIB $GRBFILE | grep "HGT"| grep "cld base" | awk -F ":" '{print $1}')
            export index_t=$($WGRIB $GRBFILE | grep "HGT"| grep "cld top" | awk -F ":" '{print $1}')
            time $GDALCALC -A $GRBFILE --A_band $index_t -B $GRBFILE --B_band $index_b --outfile $param.tiff --calc="A-B"
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "VISIB")
            export index=$($WGRIB $GRBFILE | grep "RH" | grep "2 m above" | awk -F ":" '{print $1}')
            time $GDALCALC -A $GRBFILE --A_band $index --outfile $param.tiff --calc="(A>74)*(192.130 -41.7*numpy.log(A+1.) + 0.35)" --NoDataValue=0
            time $GDALCALC -A $param.tiff --A_band 1 --outfile $param.tif --calc="A*(A<10)" --NoDataValue=0
	    time $GDALPOLYGON -b 1 -of 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "LWRAD")
            export index=$($WGRIB $GRBFILE | grep "ULWRF" | grep "sfc" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "TROP")
            export index=$($WGRIB $GRBFILE | grep "TMP" | grep "tropopause" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "FRZH")
            export index=$($WGRIB $GRBFILE | grep "HGT" | grep "0C isotherm:" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "TMP850")
            export index=$($WGRIB $GRBFILE | grep "TMP" | grep "850 mb" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "TMP500")
            export index=$($WGRIB $GRBFILE | grep "TMP" | grep "500 mb" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "GP850C")
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
	    /bin/rm  $param.shp $param.gs

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
	    /bin/rm  $param.shp $param.gs

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
        "GP850")
            export index=$($WGRIB $GRBFILE | grep "HGT" | grep "850 mb" | awk -F ":" '{print $1}')
	    time $GDALPOLYGON -b $index -of 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME
            ;;
        "DIR850")
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
        *)
            echo "Unrecognized parameter $param!"
            continue
            ;;
    esac
  fi

  if [ "$param" != "GP850C" ] && [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
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
