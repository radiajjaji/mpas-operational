#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --job-name=JSN.MPS.3
#SBATCH --ntasks=8
#SBATCH --partition=opr,workq --no-requeue
#SBATCH --time=3:00:00
#SBATCH --propagate=NONE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
. /scratch/lus/arw/model/slurm/modelenv.sh

export ATTRNAME='DN'
export WIND_THIN_FACTOR=4
export THIN_FACTOR=10

cd $TMPDIR

echo "Running on :" $(hostname)

export CYCLE=$(expr $START_DATE | cut -c9-10)
export RUN_DATE=$(expr $START_DATE | cut -c1-8)

export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)

export GRBFILE=$D/r$CYCLE/mpsprs.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.grb1
export MPSOUTFILE=$D/r$CYCLE/mpsout.d01.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.nc

[ -f $GRBFILE ] || exit 0

export parameters=("LITN" "DBZ" "BRTMP" "APCP6H" "GP700C" "GP700" "DIR700" "GP500" "GP500C" "DIR500" "SPDJET" "DIRJET" "CAPE" "CIN" "LINDEX" "VERVEL")

# --- 6H Accumulation Logic ---
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
else
    # Remove APCP6H if not on a 6-hour interval
    export parameters=("LITN" "DBZ" "BRTMP" "GP700C" "GP700" "DIR700" "GP500" "GP500C" "DIR500" "SPDJET" "DIRJET" "CAPE" "CIN" "LINDEX" "VERVEL")
fi
# --- End 6H Accumulation Logic ---

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

for param in "${parameters[@]}"; do

    filesize=0
    [ -f "$OUTPUTPATH/$param.geojson" ] && filesize=$(stat -c %s "$OUTPUTPATH/$param.geojson" 2>/dev/null || echo 0)
    if [ ! -f "$OUTPUTPATH/$param.geojson" ] || [ "$filesize" -le 111111 ] ; then
        [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

        case $param in
            "APCP6H")
		export CLEVS="0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0"
                if [ $SLURM_ARRAY_TASK_ID -ge 6 ] ; then
                    export index_current=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
                    export index_past=$($WGRIB $GRBFILE_6H | grep "APCP" | awk -F ":" '{print $1}')
                    time $GDALCALC -A $GRBFILE --A_band $index_current -B $GRBFILE_6H --B_band $index_past --outfile $param.tiff --calc="A-B"
                else
                    export index_current=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
                    $WGRIB -d $index_current -grib -o $param.grb $GRBFILE
                    $GDAL_TRANSLATE -of GTiff $param.grb $param.tiff
                    /bin/rm $param.grb
                fi
                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                /bin/rm $param.tiff
                ;;
            "DBZ")
                export CINT=2

                # REFD_MAX comes directly from mpsout NetCDF
                # Convert it first to GeoTIFF, then let raster_to_geojson work as usual
                time $GDAL_TRANSLATE \
                    -of GTiff \
                    -a_srs EPSG:4326 \
                    -a_ullr -180 90 180 -90 \
                    NETCDF:"$MPSOUTFILE":REFD_MAX \
                    $param.tiff

                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson \
                    -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -i $CINT

                /bin/rm -f $param.tiff
                ;;
            "DBZO")
                export index=$($WGRIB $GRBFILE | grep "REFC" | awk -F ":" '{print $1}')
                export CLEVS="5,10,20,30,35,40,45,50,55,60,65"
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                ;;
            "BRTMP")
                export index=$($WGRIB $GRBFILE | grep "BRTMP" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "GP700"|"GP500")
                export index=$($WGRIB $GRBFILE | grep "HGT" | grep "${param:2:3} mb" | awk -F ":" '{print $1}')
                export CINT=20
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "SPDJET")
                export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "250 mb" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "250 mb" | awk -F ":" '{print $1}')
                export CINT=5
                time $GDALCALC -A $GRBFILE --A_band $index_u -B $GRBFILE --B_band $index_v --outfile $param.tiff \
                    --calc="(1.944*sqrt(A*A+B*B))*((1.944*sqrt(A*A+B*B)) > 50)"
                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -i $CINT
                /bin/rm $param.tiff
                ;;
            "CAPE")
                export index=$($WGRIB $GRBFILE | grep "CAPE" | awk -F ":" '{print $1}')
                export CINT=200
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "CIN")
                export index=$($WGRIB $GRBFILE | grep "CIN" | awk -F ":" '{print $1}')
                export CINT=200
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "LINDEX")
                export index=$($WGRIB $GRBFILE | grep "4LFTX" | awk -F ":" '{print $1}')
                export CINT=0.5
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "VERVEL")
                export index=$($WGRIB $GRBFILE | grep "VVEL" | grep "850 mb" | awk -F ":" '{print $1}')
                export CINT=0.2
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;

            "DIR700"|"DIR500"|"DIRJET")
                PRESS="${param:3:3} mb"
                if [ "$param" = "DIRJET" ]; then PRESS="250 mb"; fi
                
                export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "$PRESS" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "$PRESS" | awk -F ":" '{print $1}')
                $WGRIB -d $index_u -grib -o u.grb $GRBFILE
                $WGRIB -d $index_v -grib -o v.grb $GRBFILE
                
                # Calculate Speed (SPD) and Direction (DIR)
                $CDO -nint -mulc,1.994 -sqrt -add -sqr u.grb -sqr v.grb spd.grb
                # DIR = 180 + 57.3 * atan2(-U, -V)
                $CDO -nint -addc,180. -mulc,57.3 -atan2 -mulc,-1 u.grb -mulc,-1 v.grb dir.grb
                
                # Combine SPD and DIR into a single TIFF: (SPD * 1000) + DIR
                time $GDALCALC -A spd.grb -B dir.grb --outfile $param.tiff --calc="A*1000 + B"

                # Raster to GeoJSON points, thinning factor applied
                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pt -b 1 -a $ATTRNAME -f 'GeoJSON' -s $WIND_THIN_FACTOR

                /bin/rm $param.tiff u.grb v.grb spd.grb dir.grb
                ;;

            "GP700C"|"GP500C")
                PRESS="${param:2:3} mb"
                export index=$($WGRIB $GRBFILE | grep "HGT" | grep "$PRESS" | awk -F ":" '{print $1}')
                export CINT=20

                # Raster to GeoJSON lines (ln) and points (pt), thinning factor applied
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t ln pt -s $THIN_FACTOR -a $ATTRNAME -f 'GeoJSON' -i $CINT -b $index 
                ;;

            "LITN")
                # Identify indexes for Reflectivity and CAPE
                export idx_ref=$($WGRIB $GRBFILE | grep "REFC" | awk -F ":" '{print $1}')
                export idx_cape=$($WGRIB $GRBFILE | grep "CAPE" | awk -F ":" '{print $1}')

                # Formula: If Reflectivity > 35 AND CAPE > 500, output Reflectivity value, else 0
                # This highlights the "core" of the storm capable of lightning
                time $GDALCALC -A $GRBFILE --A_band $idx_ref -B $GRBFILE --B_band $idx_cape \
                    --outfile $param.tiff \
                    --calc="1.8 * A * (A > 35) * (B > 500)" --NoDataValue=0

                # Define levels of intensity (mapping back to dBZ strength in the core)
                export CLEVS="35,40,45,50,55,60"

                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                /bin/rm $param.tiff
                ;;

        esac
    fi

    mbtilesize=0
    [ -f "$OUTPUTPATH/$param.mbtiles" ] && mbtilesize=$(stat -c %s "$OUTPUTPATH/$param.mbtiles" 2>/dev/null || echo 0)
    if [ ! -f "$OUTPUTPATH/$param.mbtiles" ] || [ "$mbtilesize" -le 124576 ] ; then
        [ -f $OUTPUTPATH/$param.geojson ] || continue

        [ -f "$OUTPUTPATH/$param.geojson" ] || continue
        [ -f "$OUTPUTPATH/$param.geojson" ] || continue
        time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o "$OUTPUTPATH/$param.mbtiles" "$OUTPUTPATH/$param.geojson"
        
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
