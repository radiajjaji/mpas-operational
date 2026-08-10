#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --job-name=JSN.MPS.2
#SBATCH --ntasks=8
#SBATCH --partition=opr,workq --no-requeue
#SBATCH --time=3:00:00
#SBATCH --propagate=NONE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

export SLURM_IGNORE_LIMITS=1 || true

set -x
echo $(hostname) >> /scratch/lus/arw/model/log/json.mpas.part2.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh

export WIND_THIN_FACTOR=2
export ATTRNAME='DN'
export THIN_FACTOR=10

cd $TMPDIR
set +x

module load PrgEnv-intel
module load cray-netcdf
module load cray-hdf5
module load craype-x86-rome
module load cray-libsci

ulimit -s unlimited
ulimit -c 0
set -x

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

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/
mkdir -p "$OUTPUTPATH"
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

declare -a params=("GUST" "LCLD" "MCLD" "HCLD" "CBCLD" "VISIB" "LWRAD" "TROP" "FRZH" "TMP850" "TMP500" "GP850C" "GP850" "DIR850")

for param in "${params[@]}"; do
    /bin/rm -f "$param.tiff" "$param.tif" u.grb v.grb spd.grb dir.grb

    filesize=0
    [ -f "$OUTPUTPATH/$param.geojson" ] && filesize=$(stat -c %s "$OUTPUTPATH/$param.geojson" 2>/dev/null || echo 0)

    # Check if GeoJSON needs to be created
    if [ ! -f "$OUTPUTPATH/$param.geojson" ] || [ "$filesize" -le 111111 ]; then
        [ ! -f "$OUTPUTPATH/$param.geojson" ] || /bin/rm "$OUTPUTPATH/$param.geojson"

        case $param in
            # --- Parameters Requiring Isolines (Polygons: -t pl) ---
            "GUST")
                export index=$($WGRIB "$GRBFILE" | grep "GUST" | awk -F ":" '{print $1}')
                export CINT=5
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "LCLD"|"MCLD"|"HCLD")
                if [ "$param" = "LCLD" ]; then pattern="LCDC";
                elif [ "$param" = "MCLD" ]; then pattern="MCDC";
                elif [ "$param" = "HCLD" ]; then pattern="HCDC"; fi
                export index=$($WGRIB "$GRBFILE" | grep "$pattern" | awk -F ":" '{print $1}')
                export CINT=5
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "LWRAD")
                export index=$($WGRIB "$GRBFILE" | grep "ULWRF" | grep "sfc" | awk -F ":" '{print $1}')
                export CINT=25
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "TROP")
                export index=$($WGRIB "$GRBFILE" | grep "TMP" | grep "tropopause" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "FRZH")
                export index=$($WGRIB "$GRBFILE" | grep "HGT" | grep "0C isotherm:" | awk -F ":" '{print $1}')
                export CINT=200
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "TMP850")
                export index=$($WGRIB "$GRBFILE" | grep "TMP" | grep "850 mb" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "TMP500")
                export index=$($WGRIB "$GRBFILE" | grep "TMP" | grep "500 mb" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            "GP850")
                export index=$($WGRIB "$GRBFILE" | grep "HGT" | grep "850 mb" | awk -F ":" '{print $1}')
                export CINT=20
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                ;;
            
            # --- Parameters Requiring Pre-processing (GDALCALC) ---
            "CBCLD")
                export index_b=$($WGRIB "$GRBFILE" | grep "HGT" | grep "cld base" | awk -F ":" '{print $1}')
                export index_t=$($WGRIB "$GRBFILE" | grep "HGT" | grep "cld top" | awk -F ":" '{print $1}')
                export CINT=500
                time $GDALCALC -A "$GRBFILE" --A_band "$index_t" -B "$GRBFILE" --B_band "$index_b" --outfile "$param.tiff" --calc="A-B"
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm -f "$param.tiff"
                ;;
            "VISIB")
                export index=$($WGRIB "$GRBFILE" | grep "RH" | grep "2 m above" | awk -F ":" '{print $1}')
                export CINT=1.0
                time $GDALCALC -A "$GRBFILE" --A_band "$index" --outfile "$param.tiff" --calc="(A>74)*(192.130 -41.7*numpy.log(A+1.) + 0.35)" --NoDataValue=0
                time $GDALCALC -A "$param.tiff" --A_band 1 --outfile "$param.tif" --calc="A*(A<10)" --NoDataValue=0
                time $RASTER_TO_GEOJSON "$param.tif" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm -f "$param.tiff" "$param.tif"
                ;;

            # --- Parameters Requiring Multiple Feature Types or Complex Calculations ---
            "GP850C")
                export index=$($WGRIB "$GRBFILE" | grep "HGT" | grep "850 mb" | awk -F ":" '{print $1}')
                export CINT=20
                # Extract line (ln) and thinned point (pt) features for Contours
                time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t ln pt -s $THIN_FACTOR -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT" -b "$index" -s 20
                ;;
            "DIR850")
                export index_u=$($WGRIB "$GRBFILE" | grep "UGRD" | grep "850 mb" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB "$GRBFILE" | grep "VGRD" | grep "850 mb" | awk -F ":" '{print $1}')
                
                $WGRIB -d "$index_u" -grib -o u.grb "$GRBFILE"
                $WGRIB -d "$index_v" -grib -o v.grb "$GRBFILE"

                # Calculate Speed (SPD) and Direction (DIR) using CDO
                $CDO -nint -mulc,1.994 -sqrt -add -sqr u.grb -sqr v.grb spd.grb
                $CDO -nint -addc,180. -mulc,57.3 -atan2 -mulc,-1 u.grb -mulc,-1 v.grb dir.grb
                
                # Combine SPD and DIR into a single TIFF: (SPD * 1000) + DIR
                time $GDALCALC -A spd.grb -B dir.grb --outfile "$param.tiff" --calc="A*1000 + B"

                # Extract thinned point (pt) features for Wind Barbs/Arrows
                time $RASTER_TO_GEOJSON "$param.tiff" "$OUTPUTPATH/$param.geojson" -t pt -b 1 -a "$ATTRNAME" -f 'GeoJSON' -s "$WIND_THIN_FACTOR"

                /bin/rm -f "$param.tiff" u.grb v.grb spd.grb dir.grb
                ;;
            *)
                continue
                ;;
        esac
    fi

    # MBTiles processing for ALL parameters
    mbtilesize=0
    [ -f "$OUTPUTPATH/$param.mbtiles" ] && mbtilesize=$(stat -c %s "$OUTPUTPATH/$param.mbtiles" 2>/dev/null || echo 0)
    if [ "$mbtilesize" -le 124576 ]; then
        # Check if geojson exists before running tpcanoe
        [ -f "$OUTPUTPATH/$param.geojson" ] || continue

        time $TPCANOE --layer='MPS_'$param --name='MPS_'$param -o "$OUTPUTPATH/$param.mbtiles" "$OUTPUTPATH/$param.geojson"
        CNT=$(sqlite3 "$OUTPUTPATH/$param.mbtiles" "select count(*) from tiles;" 2>/dev/null || echo 0)
        if [ "$CNT" -le 0 ]; then
          echo "SKIP $param: empty tiles table => $OUTPUTPATH/$param.mbtiles"
        else
          rsync -avz --delay-updates --partial-dir=".rsync-partial" "$OUTPUTPATH/$param.mbtiles" $USER@$WEB:$REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/
        fi
    fi
done
