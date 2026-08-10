#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --job-name=JSN.MPS.1
#SBATCH --ntasks=8
#SBATCH --partition=opr,workq --no-requeue
#SBATCH --time=3:00:00
#SBATCH --propagate=NONE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname) >>/scratch/lus/arw/model/log/json.mpas.part1.$SLURM_ARRAY_TASK_ID.out

. /scratch/lus/arw/model/slurm/modelenv.sh

export ATTRNAME='DN'
export WIND_THIN_FACTOR=2
export THIN_FACTOR=10

cd "$TMPDIR"

echo "Running on :" $(hostname)

export CYCLE=$(expr $START_DATE | cut -c9-10)
export RUN_DATE=$(expr $START_DATE | cut -c1-8)

export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE | cut -c9-10)

export GRBFILE=$D/r$CYCLE/mpsprs.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.grb1
[ -f "$GRBFILE" ] || exit 0

# --- 1H Date Calculation Logic ---
if [ $SLURM_ARRAY_TASK_ID -gt 0 ] ; then
    export DATE_1H=$($SMSDATE -1 $ACTUAL_DATE)
else
    export DATE_1H=$ACTUAL_DATE
fi

export YEAR_1H=$(expr $DATE_1H | cut -c1-4)
export MONTH_1H=$(expr $DATE_1H | cut -c5-6)
export DAY_1H=$(expr $DATE_1H | cut -c7-8)
export HOUR_1H=$(expr $DATE_1H | cut -c9-10)

export GRBFILE_1H=$D/r$CYCLE/mpsprs.$YEAR_1H-$MONTH_1H-${DAY_1H}_${HOUR_1H}_00_00.grb1
# --- End 1H Date Calculation Logic ---

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/
mkdir -p "$OUTPUTPATH"
ssh "$USER@$WEB" mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

export param_list=("MSLPC" "MSLP" "TMP2M" "RH2M" "APCP1H" "SPD10M" "DIR10M" "WAVES")

for param in "${param_list[@]}"; do

    filesize=0
    [ -f "$OUTPUTPATH/$param.geojson" ] && filesize=$(stat -c %s "$OUTPUTPATH/$param.geojson" 2>/dev/null || echo 0)

    if [ ! -f "$OUTPUTPATH/$param.geojson" ] || [ "$filesize" -le 111111 ] ; then
        [ ! -f "$OUTPUTPATH/$param.geojson" ] || /bin/rm "$OUTPUTPATH/$param.geojson"

        case $param in
            "MSLP")
                # Mean Sea Level Pressure (PRMSL) to hPa (0.01 * A)
                export index=$($WGRIB "$GRBFILE" | grep "PRMSL" | awk -F ":" '{print $1}')
                export CINT=2
                time $GDALCALC -A "$GRBFILE" --A_band="$index" --calc="0.01*A" --outfile="$param.tiff"
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm "$param.tiff"
                ;;
            "TMP2M")
                # Temperature at 2m (TMP)
                export index=$($WGRIB "$GRBFILE" | grep "TMP" | grep "2 m above gnd" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "RH2M")
                # Relative Humidity at 2m (RH) with CDO smoothing and thresholding
                export index=$($WGRIB "$GRBFILE" | grep "RH" | grep "2 m above gnd" | awk -F ":" '{print $1}')
		export CLEVS="70,72,74,76,78,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100"
                $WGRIB -d "$index" -grib -o "$param.grb" "$GRBFILE"
                $CDO -smooth9 "$param.grb" "$param.cdo.grb"
                time $GDALCALC -A "$param.cdo.grb" --A_band=1 --calc="A*(A>74)" --NoDataValue=0 --outfile="$param.tiff"
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -lev "$CLEVS"
                /bin/rm "$param.grb" "$param.cdo.grb" "$param.tiff"
                ;;
            "APCP1H")
                # 1-Hour Accumulated Precipitation (APCP) difference
                export index=$($WGRIB "$GRBFILE" | grep "APCP" | awk -F ":" '{print $1}')
                export CLEVS="0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,14.0,18.0,22.0,30.0,35.0,40.0"
                time $GDALCALC -A "$GRBFILE" --A_band="$index" -B "$GRBFILE_1H" --B_band="$index" --outfile="$param.tiff" --calc="A-B"
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -levs "$CLEVS"
                /bin/rm "$param.tiff"
                ;;
            "WAVES")
                # Calculated Significant Wave Height (WAVES) based on 10m wind speed
                export index_u=$($WGRIB "$GRBFILE" | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB "$GRBFILE" | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
                export index_mask=$($WGRIB "$GRBFILE" | grep "LAND" | awk -F ":" '{print $1}')
                export CINT=0.5
                # Calculate wind speed (A): 1.994 * sqrt(U^2 + V^2)
                time $GDALCALC -A "$GRBFILE" --A_band="$index_u" -B "$GRBFILE" --B_band="$index_v" \
                    --outfile="$param.tiff" --calc="1.994*numpy.sqrt(A*A+B*B)"
                # Calculate wave height (B): (0.071 + 0.232*A + 0.009*A*A) * (LAND < 0.5)
                time $GDALCALC -A "$param.tiff" --A_band=1 -B "$GRBFILE" --B_band="$index_mask" --outfile="$param.tif" \
                    --calc="(0.071 + 0.232*A + 0.009*A*A)*(B<0.5)" --NoDataValue=0
                time $RASTER_TO_GEOJSON "$param.tif" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm "$param.tiff" "$param.tif"
                ;;
            "SPD10M")
                # 10m Wind Speed calculation
                export index_u=$($WGRIB "$GRBFILE" | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB "$GRBFILE" | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
                export CINT=2
                time $GDALCALC -A "$GRBFILE" --A_band="$index_u" -B "$GRBFILE" --B_band="$index_v" \
                    --outfile="$param.tiff" --calc="1.994*numpy.sqrt(A*A+B*B)"
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm "$param.tiff"
                ;;
            "DIR10M")
                # 10m Wind Direction (DIR) with thinning
                export index_u=$($WGRIB "$GRBFILE" | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB "$GRBFILE" | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')

                $WGRIB -d "$index_u" -grib -o u.grb "$GRBFILE"
                $WGRIB -d "$index_v" -grib -o v.grb "$GRBFILE"

                # Calculate Speed (SPD) and Direction (DIR) using CDO
                $CDO -nint -mulc,1.994 -sqrt -add -sqr u.grb -sqr v.grb spd.grb
                $CDO -nint -addc,180. -mulc,57.3 -atan2 -mulc,-1 u.grb -mulc,-1 v.grb dir.grb

                # Combine SPD and DIR into a single TIFF: (SPD * 1000) + DIR
                time $GDALCALC -A spd.grb -B dir.grb --outfile "$param.tiff" --calc="A*1000 + B"

                # Extract thinned point (pt) features for Wind Barbs/Arrows
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pt -b 1 -a "$ATTRNAME" -f 'GeoJSON' -s "$WIND_THIN_FACTOR"

                /bin/rm "$param.tiff" u.grb v.grb spd.grb dir.grb
                ;;
            "MSLPC")
                # MSLP Contours (ln) and thinned points (pt)
                export index=$($WGRIB "$GRBFILE" | grep "PRMSL" | awk -F ":" '{print $1}')
                export CINT=1

                # MSLP to hPa
                time $GDALCALC -A "$GRBFILE" --A_band="$index" --calc="0.01*A" --outfile="$param.tiff"
                
                # Extract line (ln) and thinned point (pt) features
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t ln pt -s $THIN_FACTOR -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT" -b 1
                
                /bin/rm "$param.tiff"
                ;;
        esac
    fi

    # MBTiles processing
    mbtilesize=0
    [ -f "$OUTPUTPATH/$param.mbtiles" ] && mbtilesize=$(stat -c %s "$OUTPUTPATH/$param.mbtiles" 2>/dev/null || echo 0)
    
    if [ "$mbtilesize" -le 124576 ] ; then
        [ -f "$OUTPUTPATH/$param.geojson" ] || continue

        time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o "$OUTPUTPATH/$param.mbtiles" "$OUTPUTPATH/$param.geojson"
        CNT=$(sqlite3 "$OUTPUTPATH/$param.mbtiles" "select count(*) from tiles;" 2>/dev/null || echo 0)
        if [ "$CNT" -le 0 ]; then
          echo "SKIP $param: empty tiles table => $OUTPUTPATH/$param.mbtiles"
        else
          rsync -avz --delay-updates --partial-dir=".rsync-partial" "$OUTPUTPATH/$param.mbtiles" "$USER@$WEB:$REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/"
        fi
    fi
done
