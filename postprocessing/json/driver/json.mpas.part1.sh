#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part1.%a.out
#SBATCH --job-name=JSON.M1
#SBATCH --ntasks=8
#SBATCH --partition=uan,hpclm,slm --oversubscribe  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part1.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh


export WIND_THIN_FACTOR=2
export POLYGON_TYPE=BANDS
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

if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $((SLURM_ARRAY_TASK_ID % 24)) -eq 0 ] ; then
   export DATE_24H=$($SMSDATE -24 $ACTUAL_DATE)
   export DO_24H="yes"
else
   export DO_24H="no"
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

if [ ! -f $D/r$CYCLE/$(basename $GRBFILE_1H .grb1).ctl ] ; then
   $GRIB2CTL -verf $GRBFILE_1H > mpsprsp.ctl
   $GRIBMAP -big -i mpsprsp.ctl
else
   cp $D/r$CYCLE/$(basename $GRBFILE_1H .grb1).ctl mpsprsp.ctl
fi

# Input And Output Directories

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$START_DATE/$ACTUAL_DATE/

if [ ! -d $OUTPUTPATH ] ; then
     mkdir -p $OUTPUTPATH
fi

# Json and Vectoriel Tiles Creation for a set of selected parameters

export param='MSLPC'
if [ ! -f $OUTPUTPATH/$param.Line.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.Line.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.Line.geojson ] || /bin/rm $OUTPUTPATH/$param.Line.geojson $OUTPUTPATH/$param.Point.geojson
   export index=$($WGRIB $GRBFILE | grep "PRMSL" | awk -F ":" '{print $1}')
   #$WGRIB -d $index -grib -o mslp.grb $GRBFILE
   #$CDO -P 8 smooth,maxpoints=100 mslp.grb mslp_smooth.grb ; rm mslp.grb
   #$CDO -P 8 smooth,maxpoints=100 mslp_smooth.grb mslp.grb ; rm mslp_smooth.grb
   #export CLEVS="950 955 960 965 970 975 980 982 984 986 988 990 992 994 996 998 1000 1002 1004 1006 1008 1010 1012 1014 1016 1018 1020 1022 1024 1026 1028 1030 1032 1034 1036 1038 1040"
   #time $GDALCALC -A mslp.grb --A_band=$index --calc="0.01*A" --outfile=$param.tiff ; rm mslp.grb
   #time $GDALCONT $param.tiff   -b 1 -a $ATTRNAME -f 'GeoJSON' -fl $CLEVS -q $OUTPUTPATH/${param}.Line.geojson
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set cint 2'
'd smth9(smth9(smth9(prmslm --oversubscribesl/100.)))'
'quit'
EOF

time $GRADS -lbc '$param.gs' ; /bin/rm $param.gs
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Line.geojson $param.shp ; sed -i s"/CNTR_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Line.geojson

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -pt $param -fmt 6 1'
'define pressure=smth9(smth9(smth9(prmslm --oversubscribesl/100.)))'
'set cint 2'
'd skip(pressure,20)'
'quit'
EOF

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Point.geojson $param.shp ; sed -i s"/GRID_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Point.geojson

fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Line.geojson
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.Point.mbtiles $OUTPUTPATH/$param.Point.geojson
   time $TILEJOIN -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.Line.mbtiles $OUTPUTPATH/$param.Point.mbtiles
fi

export param='MSLP'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "PRMSL" | awk -F ":" '{print $1}')
   export CINT=2
   time $GDALCALC -A $GRBFILE --A_band=$index --calc="0.01*A" --outfile=$param.tiff
   time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT  $OUTPUTPATH/$param.geojson
fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='TMP2M'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "TMP" | grep "2 m above gnd" |  awk -F ":" '{print $1}')
   export CINT=2
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='RH2M'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "RH" | grep "2 m above gnd" |  awk -F ":" '{print $1}')
   if [ $POLYGON_TYPE = 'BANDS' ] ; then
      export CINT=2
      $WGRIB -d $index -grib -o $param.grb $GRBFILE
      $CDO -smooth9 $param.grb $param.cdo.grb
      time $GDALCALC -A $param.cdo.grb --A_band=1 --calc="A*(A>74)" --NoDataValue=0 --outfile=$param.tiff
      time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALCALC -A $GRBFILE --A_band=$index --calc="A*(A>74)" --NoDataValue=0 --outfile=$param.tiff
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='APCP1H'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
   if [ $POLYGON_TYPE = 'BANDS' ] ; then

cat << EOF > $param.gs
'open mpsprsp.ctl'
'define rrp=apcpsfc'
'close 1'
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set clevs 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0 3.0 4.0 5.0 6.0 7.0 8.0 10.0 15.0 20.0 30.0'
'd apcpsfc - rrp'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      #export CINT=0.5
      #time $GDALCALC -A $GRBFILE --A_band $index -B $GRBFILE_1H --B_band $index --outfile $param.tiff --calc="A-B"
      #time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='SPD10M'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
   export CINT=2

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 2'
'd 1.994*mag(ugrd10m,vgrd10m)'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

   #time $GDALCALC -A $GRBFILE --A_band $index_u -B $GRBFILE --B_band $index_v --outfile $param.tiff --calc="1.994*sqrt(A*A+B*B)"
   #time $ISOBANDS $param.tiff  -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DIR10M'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

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

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson

fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='WAVES'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "10 m above gnd" | awk -F ":" '{print $1}')
   export index_mask=$($WGRIB $GRBFILE | grep "LAND" | awk -F ":" '{print $1}')
   time $GDALCALC -A $GRBFILE --A_band $index_u -B $GRBFILE --B_band $index_v --outfile $param.tiff --calc="1.994*sqrt(A*A+B*B)"
   time $GDALCALC -A $param.tiff --A_band 1 -B $GRBFILE --B_band $index_mask --outfile $param.tif \
        --calc="(0.071 + 0.232*A + 0.009*A*A)*(B<0.5)" --NoDataValue=0

   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tif $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set clevs 0.01 0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0 10.0 12.0 14.0 18.0 20.0 22.0 24.0'
'define spd=1.994*mag(ugrd10m,vgrd10m)'
'd (0.071 + 0.232 * spd + 0.009 * spd * spd) * (1.0-landsfc)  '
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

      #export CINT=2
      #time $ISOBANDS $param.tiff  -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson 
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi
