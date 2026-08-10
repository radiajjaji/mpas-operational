#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --job-name=JSN.MPS.1
#SBATCH --ntasks=8
#SBATCH --partition=opr  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part1.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

export ATTRNAME='DN'
export WIND_THIN_FACTOR=2

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

if [ $SLURM_ARRAY_TASK_ID -gt 0 ] ; then
   export DATE_1H=$($SMSDATE -1 $ACTUAL_DATE)
else
   export DATE_1H=$ACTUAL_DATE
fi

export YEAR_1H=$(expr $DATE_1H  | cut -c1-4)
export MONTH_1H=$(expr $DATE_1H | cut -c5-6)
export DAY_1H=$(expr $DATE_1H   | cut -c7-8)
export HOUR_1H=$(expr $DATE_1H  | cut -c9-10)
export GRBFILE_1H=$D/r$CYCLE/mpsprs.$YEAR_1H-$MONTH_1H-${DAY_1H}_${HOUR_1H}_00_00.grb1

$GRIB2CTL -verf $GRBFILE_1H > mpsprsp.ctl
$GRIBMAP -big -i mpsprsp.ctl

# Input And Output Directories

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

export param_list=("MSLPC" "MSLP" "TMP2M" "RH2M" "APCP1H" "SPD10M" "DIR10M" "WAVES")

for param in "${param_list[@]}"; do

  if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
     [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

    case $param in
    "MSLP")
        export index=$($WGRIB $GRBFILE | grep "PRMSL" | awk -F ":" '{print $1}')
        time $GDALCALC -A $GRBFILE --A_band=$index --calc="0.01*A" --outfile=$param.tiff
	time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;
    "TMP2M")
        export index=$($WGRIB $GRBFILE | grep "TMP" | grep "2 m above gnd" | awk -F ":" '{print $1}')
	time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;
    "RH2M")
        export index=$($WGRIB $GRBFILE | grep "RH" | grep "2 m above gnd" | awk -F ":" '{print $1}')
	time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;
    "APCP1H")
        export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
        time $GDALCALC -A $GRBFILE --A_band=$index -B $GRBFILE_1H --B_band=$index --outfile=$param.tiff --calc="A-B"
	time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;
    "WAVES")
        export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
        export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
        export index_mask=$($WGRIB $GRBFILE | grep "LAND" | awk -F ":" '{print $1}')

        time $GDALCALC -A $GRBFILE --A_band=$index_u -B $GRBFILE --B_band=$index_v \
            --outfile=$param.tiff --calc="1.994*sqrt(A*A+B*B)"
        time $GDALCALC -A $param.tiff --A_band=1 -B $GRBFILE --B_band=$index_mask --outfile=$param.tif \
            --calc="(0.071 + 0.232*A + 0.009*A*A)*(B<0.5)" --NoDataValue=0
	time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;

    "SPD10M")
        export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
        export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
        time $GDALCALC -A $GRBFILE --A_band=$index_u -B $GRBFILE --B_band=$index_v \
            --outfile=$param.tiff --calc="1.994*sqrt(A*A+B*B)"
	time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
        ;;
    "DIR10M")
        export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
        export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')

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
'd skip(1000*ugrd10m.1+ugrd10m.2,$WIND_THIN_FACTOR)'
'quit'
EOF

	time $GRADS -lbc "$param.gs"
	time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson
        ;;
    "MSLPC")
        # Handle MSLPC
        if [ ! -f $OUTPUTPATH/$param.Line.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.Line.geojson) -le 111111 ]; then
            [ ! -f $OUTPUTPATH/$param.Line.geojson ] || /bin/rm $OUTPUTPATH/$param.Line.geojson $OUTPUTPATH/$param.Point.geojson
            export index=$($WGRIB $GRBFILE | grep "PRMSL" | awk -F ":" '{print $1}')
            
            # Create line GeoJSON
            cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set cint 1'
'd msletmsl/100.'
'quit'
EOF
            time $GRADS -lbc "$param.gs" ; /bin/rm $param.gs
            time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Line.geojson $param.shp
            sed -i s"/CNTR_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Line.geojson

            # Create point GeoJSON
            cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -pt $param -fmt 6 1'
'd skip(msletmsl/100,20)'
'quit'
EOF
            time $GRADS -lbc "$param.gs" ; /bin/rm $param.gs
            time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Point.geojson $param.shp
            sed -i s"/GRID_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Point.geojson
        fi

        # Create MBTiles
        if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ]; then
            time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Line.geojson
            time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o $OUTPUTPATH/$param.Point.mbtiles $OUTPUTPATH/$param.Point.geojson
            time $TILEJOIN --force -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Point.mbtiles
	    
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
  if [ $param != "MSLPC" ] ; then
     if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ]; then
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
  fi
done
