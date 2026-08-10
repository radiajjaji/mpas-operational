#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part3.%a.out
#SBATCH --job-name=JSON.M3
#SBATCH --ntasks=8
#SBATCH --partition=uan,hpclm,slm --oversubscribe  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1

set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part3.$SLURM_ARRAY_TASK_ID.out
. /scratch/lus/arw/model/slurm/modelenv.sh


export POLYGON_TYPE=BANDS
export WIND_THIN_FACTOR=4
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

if [ ! -f $D/r$CYCLE/$(basename $GRBFILE_6H .grb1).ctl ] ; then
   $GRIB2CTL -verf $GRBFILE_6H > mpsprsp.ctl
   $GRIBMAP -big -i mpsprsp.ctl
else
   cp $D/r$CYCLE/$(basename $GRBFILE_6H .grb1).ctl mpsprsp.ctl
fi

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

export param='DBZ'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
      if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cmin 5'
'set clevs 5 10 20 30 35 40 45 50 55 60 65'
'd refcclm'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
   fi

   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      export index=$($WGRIB $GRBFILE | grep "REFC" | awk -F ":" '{print $1}')
      time $GDALCALC -A $GRBFILE --A_band=$index --calc="A*(A>5)" --NoDataValue=0 --outfile=$param.tiff
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='BRTMP'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "BRTMP" | awk -F ":" '{print $1}')
   time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

if [ $DO_6H = 'yes' ] ; then
   export param='APCP6H'
   if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
      [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
      export index=$($WGRIB $GRBFILE | grep "APCP" | awk -F ":" '{print $1}')
      export CINT=2.

      if [ $POLYGON_TYPE = 'BANDS' ] ; then
         if [ $SLURM_ARRAY_TASK_ID -ge 6 ] ; then
cat << EOF > $param.gs
'open mpsprsp.ctl'
'define rrp=apcpsfc'
'close 1'
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set clevs 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 15.0 20.0 25.0 30.0'
'd apcpsfc - rrp'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

            #time $GDALCALC -A $GRBFILE --A_band $index -B $GRBFILE_6H --B_band $index --outfile $param.tiff --calc="A-B"
            #time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
         else
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set clevs 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0 3.0 4.0 5.0 6.0 7.0 8.0 10.0 15.0 20.0 30.0 40.0 50.0 60.0'
'd apcpsfc '
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
            #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
         fi
      fi
      if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
         if [ $SLURM_ARRAY_TASK_ID -ge 6 ] ; then
            time $GDALCALC -A $GRBFILE --A_band $index -B $GRBFILE_6H --B_band $index --outfile $param.tiff --calc="A-B"
            time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
         else
            time $GDALPOLYGON $GRBFILE -b $index -f 'GeoJSON'  $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
         fi
      fi
   fi

   if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
      time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
   fi
fi

export param='GP700C'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT" | grep "700 mb" | awk -F ":" '{print $1}')

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set lev 700'
'set cint 20'
'd smth9(smth9(smth9(hgtprs)))'
'quit'
EOF

time $GRADS -lbc '$param.gs' ; /bin/rm $param.gs
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Line.geojson $param.shp ; sed -i s"/CNTR_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Line.geojson

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -pt $param -fmt 6 1'
'set lev 700'
'define geop=smth9(smth9(smth9(hgtprs/100.)))'
'set cint 20'
'd skip(hgtprs,20)'
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

export param='GP700'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT" | grep "700 mb" | awk -F ":" '{print $1}')
   export CINT=20
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DIR700'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "700 mb" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "700 mb" | awk -F ":" '{print $1}')

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
'd skip(1000*ugrd700mb.1+ugrd700mb.2,$WIND_THIN_FACTOR)'
'quit'
EOF

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson

fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi


export param='GP500'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT" | grep "500 mb" | awk -F ":" '{print $1}')

   export CINT=20
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT  $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='GP500C'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT" | grep "500 mb" | awk -F ":" '{print $1}')

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set lev 500'
'set cint 20'
'd smth9(smth9(smth9(hgtprs)))'
'quit'
EOF

time $GRADS -lbc '$param.gs' ; /bin/rm $param.gs
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.Line.geojson $param.shp ; sed -i s"/CNTR_VALUE/$ATTRNAME/"g $OUTPUTPATH/$param.Line.geojson

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -pt $param -fmt 6 1'
'set lev 500'
'define geop=smth9(smth9(smth9(hgtprs/100.)))'
'set cint 20'
'd skip(hgtprs,20)'
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

export param='DIR500'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "500 mb" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "500 mb" | awk -F ":" '{print $1}')

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
'd skip(1000*ugrd500mb.1+ugrd500mb.2,$WIND_THIN_FACTOR)'
'quit'
EOF

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson

fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='SPDJET'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

   export index_u=$($WGRIB $GRBFILE | grep "UGRD" | grep "250 mb" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD" | grep "250 mb" | awk -F ":" '{print $1}')
   export CINT=5

   time $GDALCALC -A $GRBFILE --A_band $index_u -B $GRBFILE --B_band $index_v --outfile $param.tiff \
        --calc="(1.994*sqrt(A*A+B*B))*((1.994*sqrt(A*A+B)) > 50)"
   time $ISOBANDS $param.tiff  -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DIRJET'

if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

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

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson

fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='CAPE'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "CAPE" | awk -F ":" '{print $1}')
   export CINT=200
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='CIN'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "CIN" | awk -F ":" '{print $1}')
   export CINT=200
   time $GDALCONT $GRBFILE -b $index -a $ATTRNAME -f 'GeoJSON' -i $CINT -q $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='LINDEX'

if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "LFTX" | grep "500-1000 mb" | awk -F ":" '{print $1}')
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cmax 0'
'set cint 0.5'
'd pli30_0mb'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

   #export CINT=2
   #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='VERVEL'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "VVEL" | grep "850 mb" | awk -F ":" '{print $1}')

cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 0.2'
'd VVEL850mb'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

   #export CINT=0.5
   #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi
