#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --job-name=JSN.MPS.4
#SBATCH --ntasks=8
#SBATCH --partition=opr,workq --no-requeue
#SBATCH --time=3:00:00
#SBATCH --propagate=NONE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

export SLURM_IGNORE_LIMITS=1 || true
export THIN_FACTOR=10

set -x
echo $(hostname) >>/scratch/lus/arw/model/log/json.mpas.part4.$SLURM_ARRAY_TASK_ID.out
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
set -x

echo "Running on :" $(hostname)

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

# --- Conditional 24H date calculation logic ---
export DO_24H="no" # Initialize DO_24H

if [ "$START_HOUR" = "00" ] ; then
    if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 24)) -eq 0 ] ; then
        export DO_24H="yes"
    fi
fi
if [ "$START_HOUR" = "12" ] ; then
    if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 24)) -eq 12 ] ; then
        export DO_24H="yes"
    fi
fi

# Define parameters based on the 24H flag and calculate INC
if [ "$DO_24H" = 'yes' ] ; then
    if [ "$START_HOUR" = "00" ] ; then
        INC=24
    elif [ "$START_HOUR" = "12" -a $SLURM_ARRAY_TASK_ID -eq 12 ] ; then
        INC=12 # First 24-hr period starts at 12H for H=12 run
    else
        INC=24
    fi

    export DATE_24H=$($SMSDATE -$INC $ACTUAL_DATE)
    export YEAR_24H=$(expr $DATE_24H | cut -c1-4)
    export MONTH_24H=$(expr $DATE_24H | cut -c5-6)
    export DAY_24H=$(expr $DATE_24H | cut -c7-8)
    export HOUR_24H=$(expr $DATE_24H | cut -c9-10)
    export GRBFILE_24H=$D/r$CYCLE/mpsprs.$YEAR_24H-$MONTH_24H-${DAY_24H}_${HOUR_24H}_00_00.grb1
    
    if [ -f $GRBFILE_24H ] ; then
        export params=("APCP24H" "RH700" "OLR" "PWATER" "VOR500" "CEIL" "DUST" "DVIS")
    else
        export params=("RH700" "OLR" "PWATER" "VOR500" "CEIL" "DUST" "DVIS")
    fi
else
    export params=("RH700" "OLR" "PWATER" "VOR500" "CEIL" "DUST" "DVIS")
fi

# Input And Output Directories
export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

mkdir -p $OUTPUTPATH
ssh $USER@$WEB mkdir -p $REMOTE_GRAPH_MPAS/mbtiles/$RUN_DATE/$CYCLE/$ACTUAL_DATE/

# Process each parameter
for param in "${params[@]}"; do
    /bin/rm -f $param.tiff $param.tif $param.grb
    filesize=0
    [ -f "$OUTPUTPATH/$param.geojson" ] && filesize=$(stat -c %s "$OUTPUTPATH/$param.geojson" 2>/dev/null || echo 0)
    if [ ! -f "$OUTPUTPATH/$param.geojson" ] || [ "$filesize" -le 111111 ] ; then
        [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

        case $param in
            "APCP24H")
                export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
		export CLEVS="0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,6.0,8.0,10.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,70.0,80.0,90.0"
                if [ "$DO_24H" = 'yes' -a $SLURM_ARRAY_TASK_ID -ge $INC ] ; then
                    time $GDALCALC -A $GRBFILE --A_band=$index -B $GRBFILE_24H --B_band=$index --outfile=$param.tiff --calc="A-B"
                    time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                else
                    time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                fi
                /bin/rm -f $param.tiff
                ;;
            "RH700")
                export index=$($WGRIB $GRBFILE | grep "RH" | grep "700 mb" | awk -F ":" '{print $1}')
                export CLEVS="60,64,68,72,74,76,78,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0"
                time $GDALCALC -A $GRBFILE --A_band=$index --outfile=$param.tiff --calc="A*(A>=50)" --NoDataValue=0
                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                /bin/rm -f $param.tiff
                ;;
            "OLR")
                export index=$($WGRIB $GRBFILE | grep "ULWRF" | grep "top" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "PWATER")
                export index=$($WGRIB $GRBFILE | grep "PWAT" | awk -F ":" '{print $1}')
                export CINT=2
                time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "CEIL")
                export index=$($WGRIB $GRBFILE | grep "HGT" | grep "cld base" | awk -F ":" '{print $1}')
                export CINT=2
                time $GDALCALC -A $GRBFILE --A_band=$index --calc="(A/1000.)" --outfile=$param.tiff
                time $RASTER_TO_GEOJSON $param.tiff $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b 1 -a $ATTRNAME -f 'GeoJSON' -i $CINT
                /bin/rm -f $param.tiff
                ;;
	    "VOR500")
                export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "500 mb" | awk -F ":" '{print $1}')
                export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "500 mb" | awk -F ":" '{print $1}')
                export CINT=2
          
                # 1. Extract U and V components into temporary GRIB files
                $WGRIB -d "$index_u" -grib -o u.grb "$GRBFILE"
                $WGRIB -d "$index_v" -grib -o v.grb "$GRBFILE"
                
                # 2. Merge U and V into one file for uv2dv (CDO expects both in one file)
                $CDO merge u.grb v.grb uv.grb
         
                # 3. Calculate Divergence (div) and Vorticity (vor), then select only Vorticity
                # Vorticity is often scaled by 10^5 or 10^4 (as you had) for typical values
                # Using -setname is crucial for CDO to correctly identify the output as 'vor'
                $CDO -setname,vor -nint -mulc,10000 -selvar,vor -uv2dv uv.grb "$param.grb"
        
                time $RASTER_TO_GEOJSON "$param.grb" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b 1 -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                /bin/rm -f u.grb v.grb uv.grb "$param.grb"
                ;;
            "DUST")
                   export index=$($WGRIB $GRBFILE | grep "PPNN" | awk -F ":" '{print $1}')
                   export CINT=500.
                   time $RASTER_TO_GEOJSON --smooth $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b ${index} -a $ATTRNAME -f 'GeoJSON' -i $CINT
                   ;;

             "DVIS")
                   export index=$($WGRIB $GRBFILE | grep "VIS" |grep "top" | awk -F ":" '{print $1}')
                   export CINT=100
                   time $RASTER_TO_GEOJSON "$GRBFILE" "$OUTPUTPATH/$param.geojson" -t pl pt -s $THIN_FACTOR -b "$index" -a "$ATTRNAME" -f 'GeoJSON' -i "$CINT"
                   ;;
            "DUST1")
                export index=$($WGRIB $GRBFILE | grep "PPNN" | awk -F ":" '{print $1}')
                export CINT=500.
                time $RASTER_TO_GEOJSON --smooth $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT
                ;;
            "DVIS1")
                export index=$($WGRIB $GRBFILE | grep "PPNN" | awk -F ":" '{print $1}')
                gdal_calc.py -A $GRBFILE --A_band=$index \
			--calc="minimum(10000, maximum(100, 10000.0/(1.0 + 0.002*((A/15.)**0.8))))" \
                  --outfile=DVIS.tif \
                  --type=Float32 \
                  --NoDataValue=10000 \
                  --quiet
                export CLEVS="10000,9000,8000,7000,6000,5000,4000,3000,2500,2000,1500,1000,800,600,500,400,300,200,100"
                time $RASTER_TO_GEOJSON --smooth $param.tif $OUTPUTPATH/$param.geojson  -t pl pt -b 1  -a $ATTRNAME -s $THIN_FACTOR  -f 'GeoJSON'  -levs $CLEVS
		/bin/rm DVIS.tif
                ;;

                #export index=$($WGRIB $GRBFILE | grep "VIS" | grep cld | awk -F ":" '{print $1}')
                #export CLEVS="10000,9000,8000,7000,6000,5000,4000,3000,2500,2000,1500,1000,800,600,500,400,300,200,100,50"
                #time $RASTER_TO_GEOJSON $GRBFILE $OUTPUTPATH/$param.geojson -t pl pt -s $THIN_FACTOR -b $index -a $ATTRNAME -f 'GeoJSON' -levs $CLEVS
                #;;

        esac
    fi

    mbtilesize=0
    [ -f "$OUTPUTPATH/$param.mbtiles" ] && mbtilesize=$(stat -c %s "$OUTPUTPATH/$param.mbtiles" 2>/dev/null || echo 0)
    if [ ! -f "$OUTPUTPATH/$param.mbtiles" ] || [ "$mbtilesize" -le 124576 ] ; then
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
