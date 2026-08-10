#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part4.%a.out
#SBATCH --job-name=JSON.M4
#SBATCH --ntasks=8
#SBATCH --partition=uan,hpclm,slm --oversubscribe  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part4.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh


export POLYGON_TYPE=BANDS
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
[ ! -f $MODEL_RUN/r$CYCLE/OK.JSON.MPAS.$START_DATE ] || \
/bin/rm $MODEL_RUN/r$CYCLE/OK.JSON.MPAS.$START_DATE

# Calculate the different dates and times

export START_YEAR=$(expr $START_DATE  | cut -c1-4)
export START_MONTH=$(expr $START_DATE | cut -c5-6)
export START_DAY=$(expr $START_DATE   | cut -c7-8)
export START_HOUR=$(expr $START_DATE  | cut -c9-10)

export ACTUAL_DATE=$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)
export ACTUAL_YEAR=$(expr $ACTUAL_DATE  | cut -c1-4)
export ACTUAL_MONTH=$(expr $ACTUAL_DATE | cut -c5-6)
export ACTUAL_DAY=$(expr $ACTUAL_DATE   | cut -c7-8)
export ACTUAL_HOUR=$(expr $ACTUAL_DATE  | cut -c9-10)
export GRBFILE=$D/r$CYCLE/mpsprs.$ACTUAL_YEAR-$ACTUAL_MONTH-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.grb1

[ -f $GRBFILE ] || exit 0

if [ ! -f $D/r$CYCLE/$(basename $GRBFILE .grb1).ctl ] ; then
   $GRIB2CTL -verf $GRBFILE > mpsprs.ctl
   $GRIBMAP -big -i mpsprs.ctl
else
   cp $D/r$CYCLE/$(basename $GRBFILE .grb1).ctl mpsprs.ctl
fi

if [ $SLURM_ARRAY_TASK_ID -gt 0 ] ; then
   export DATE_1H=$($SMSDATE -1 $ACTUAL_DATE) 
else
   export DATE_1H=$ACTUAL_DATE
fi

if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 6)) -eq 0 ] ; then
   export DATE_6H=$($SMSDATE -6 $ACTUAL_DATE)
   export DO_6H="yes"
else
   export DO_6H="no"
fi

if [ $START_HOUR = "00" ] ; then
if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 24)) -eq 0  ] ; then
   export DATE_24H=$($SMSDATE -24 $ACTUAL_DATE)
   export DO_24H="yes"
else
   export DO_24H="no"
fi
fi
if [ $START_HOUR = "12" ] ; then
if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 24)) -eq 12  ] ; then
   if [ $SLURM_ARRAY_TASK_ID -eq 12 ] ; then
      export DATE_24H=$($SMSDATE -12 $ACTUAL_DATE)
   else
      export DATE_24H=$($SMSDATE -24 $ACTUAL_DATE)
   fi
   
   export DO_24H="yes"
else
   export DO_24H="no"
fi
fi

export YEAR_1H=$(expr $DATE_1H  | cut -c1-4)
export MONTH_1H=$(expr $DATE_1H | cut -c5-6)
export DAY_1H=$(expr $DATE_1H   | cut -c7-8)
export HOUR_1H=$(expr $DATE_1H  | cut -c9-10)
export GRBFILE_1H=$D/r$CYCLE/mpsprs.$YEAR_1H-$MONTH_1H-${DAY_1H}_${HOUR_1H}_00_00.grb1

export DATE_6H=$($SMSDATE -6 $ACTUAL_DATE)
export YEAR_6H=$(expr $DATE_6H  | cut -c1-4)
export MONTH_6H=$(expr $DATE_6H | cut -c5-6)
export DAY_6H=$(expr $DATE_6H   | cut -c7-8)
export HOUR_6H=$(expr $DATE_6H  | cut -c9-10)
export GRBFILE_6H=$D/r$CYCLE/mpsprs.$YEAR_6H-$MONTH_6H-${DAY_6H}_${HOUR_6H}_00_00.grb1

export DATE_24H=$($SMSDATE -24 $ACTUAL_DATE)
export YEAR_24H=$(expr $DATE_24H  | cut -c1-4)
export MONTH_24H=$(expr $DATE_24H | cut -c5-6)
export DAY_24H=$(expr $DATE_24H   | cut -c7-8)
export HOUR_24H=$(expr $DATE_24H  | cut -c9-10)
export GRBFILE_24H=$D/r$CYCLE/mpsprs.$YEAR_24H-$MONTH_24H-${DAY_24H}_${HOUR_24H}_00_00.grb1

# Input And Output Directories

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$START_DATE/$ACTUAL_DATE/

if [ ! -d $OUTPUTPATH ] ; then
     mkdir -p $OUTPUTPATH
fi

# Json and Vectoriel Tiles Creation for a set of selected parameters

if [ $DO_24H = 'yes' ] ; then
   export param='APCP24H'
   if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
      [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
      export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
      if [ $POLYGON_TYPE = 'BANDS' ] ; then
         export CINT=1.0
         export CLEVS=" 0.2 0.4 0.6 0.8 1.0 2.0 3.0 4.0 6.0 8.0 10.0 12.0 15.0 18.0 21.0 24.0 27.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 90.0"
         if [ $SLURM_ARRAY_TASK_ID -ge 24 ] ; then
            time $GDALCALC -A $GRBFILE --A_band $index -B $GRBFILE_24H --B_band $index --outfile $param.tiff --calc="A-B"
            time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
            #time $GDALBAND $param.tiff   -b 1 -amin $ATTRNAME  -f 'GeoJSON' -fl $CLEVS -q $OUTPUTPATH/$param.geojson
         else
            time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
            #time $GDALBAND $GRBFILE   -b $index -amin $ATTRNAME  -f 'GeoJSON' -fl $CLEVS -q $OUTPUTPATH/$param.geojson
         fi
      fi
      if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
         if [ $SLURM_ARRAY_TASK_ID -ge 24 ] ; then
            time $GDALCALC -A $GRBFILE --A_band $index -B $GRBFILE_24H --B_band $index --outfile $param.tiff --calc="A-B"
            time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
         else
            time $GDALPOLYGON $GRBFILE -b $index -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
         fi
      fi
   fi
   if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
      time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
   fi
fi

export param='RH700'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "RH" | grep "700 mb" | awk -F ":" '{print $1}')
   #time $GDALCALC -A $GRBFILE --A_band=$index --calc="A*(A>60)" --NoDataValue=0 --outfile=$param.tiff

   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set lev 700'
'set cmin 50'
'set cint 1'
'd rhprs'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      #export CINT=2
      #time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='OLR'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "ULWRF" | grep "top" | awk -F ":" '{print $1}')
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 2'
'd ULWRFtoa'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

   #time $GDALPOLYGON $GRBFILE -b $index -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='PWATER'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "PWAT" | awk -F ":" '{print $1}')
   export CINT=2
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='VOR500'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "500 mb" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "500 mb" | awk -F ":" '{print $1}')

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6  $param '
'set cint 2'
'set lev 500'
'd 10000*hcurl(ugrdprs,vgrdprs)'
'quit'
EOF
   time $GRADS -lbc $param.gs
   time $OGR2OGR -f GeoJSON $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='CEIL'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT" | grep "cld base" | awk -F ":" '{print $1}')
   export CINT=2
   time $GDALCALC -A $GRBFILE --A_band=$index --calc="(A/1000.)"  --outfile=$param.tiff
   time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DUST'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "var133" | awk -F ":" '{print $1}')
   export CINT=0.1
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DVIS'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "var133" | awk -F ":" '{print $1}')
   export CINT=0.1
   time $GDALCALC -A $GRBFILE --A_band=$index --calc="12.4 * numpy.exp(-0.92 * A)"  --outfile=$param.tiff
   time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi
