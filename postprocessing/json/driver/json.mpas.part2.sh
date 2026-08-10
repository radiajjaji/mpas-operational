#!/bin/bash

#SBATCH --error=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/json.mpas.part2.%a.out
#SBATCH --job-name=JSON.M2
#SBATCH --ntasks=8
#SBATCH --partition=uan,hpclm,slm --oversubscribe  --no-requeue
#SBATCH --time=1:15:00
#SBATCH --propagate=ALL
#SBATCH --export=ALL
#SBATCH --array=0-228:1


set -x
echo $(hostname)

>>/scratch/lus/arw/model/log/json.mpas.part2.$SLURM_ARRAY_TASK_ID.out
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


# Input And Output Directories

export OUTPUTPATH=$GRAPH_MPAS/mbtiles/$START_DATE/$ACTUAL_DATE/

if [ ! -d $OUTPUTPATH ] ; then
     mkdir -p $OUTPUTPATH
fi

# Json and Vectoriel Tiles Creation for a set of selected parameters

export param='GUST'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "GUST" | awk -F ":" '{print $1}')
   export CINT=5
   time $ISOBANDS $GRBFILE   -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='LCLD'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "LCDC" | awk -F ":" '{print $1}')

   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 2'
'd lcdclcl'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      #export CINT=3
      #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='MCLD'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "MCDC" | awk -F ":" '{print $1}')
   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 2'
'd mcdcmcl'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      #export CINT=3
      #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='HCLD'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HCDC" | awk -F ":" '{print $1}')

   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set cint 2'
'd hcdchcl'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson
      #export CINT=3
      #time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b $index -f 'GeoJSON' $GRBFILE $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='CBCLD'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index_b=$($WGRIB $GRBFILE | grep "HGT"| grep "cld base" | awk -F ":" '{print $1}')
   export index_t=$($WGRIB $GRBFILE | grep "HGT"| grep "cld top" | awk -F ":" '{print $1}')
   export CINT=500
   time $GDALCALC -A $GRBFILE --A_band $index_t -B $GRBFILE --B_band $index_b --outfile $param.tiff --calc="A-B"
   time $ISOBANDS $param.tiff -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi


export param='VISIB'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "RH"| grep "2 m above" | awk -F ":" '{print $1}')
   time $GDALCALC -A  $GRBFILE --A_band $index --outfile $param.tiff --calc="(A>74)*(192.130 -41.7*numpy.log(A+1.) + 0.35)" --NoDataValue=0
   time $GDALCALC -A $param.tiff --A_band 1 --outfile $param.tif --calc="A*(A<10)" --NoDataValue=0
   if [ $POLYGON_TYPE = 'BANDS' ] ; then
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -poly -fmt 12 6 $param'
'set clevs 0.05 0.1 0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0 10.0'
*'d  60.*exp(-2.5*(rh2m-15)/80.) '
'd -0.0177*rh2m*rh2m + 1.462*rh2m + 30.8'
'quit'
EOF
      time $GRADS -lbc $param.gs
      time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/MAX_VALUE/DN/"g $OUTPUTPATH/$param.geojson

      #export CINT=0.2
      #time $ISOBANDS $param.tif -b 1 -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
   fi
   if [ $POLYGON_TYPE = 'POLYGONS' ] ; then
      time $GDALPOLYGON -b 1 -f 'GeoJSON' $param.tiff $OUTPUTPATH/$param.geojson $param $ATTRNAME -q
   fi
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='LWRAD'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "ULWRF" | grep "sfc" | awk -F ":" '{print $1}')
   export CINT=25
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT  $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='TROP'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "TMP" | grep "tropopause" | awk -F ":" '{print $1}')
   export CINT=2
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='FRZH'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT"| grep "0C isotherm:" | awk -F ":" '{print $1}')
   export CINT=200
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='TMP850'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "TMP"| grep "850 mb" | awk -F ":" '{print $1}')
   export CINT=2
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT  $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='TMP500'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "TMP"| grep "500 mb" | awk -F ":" '{print $1}')
   export CINT=2
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='GP850C'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT"| grep "850 mb" | awk -F ":" '{print $1}')
cat << EOF > $param.gs
'open mpsprs.ctl'
'set gxout shp'
'set shp -ln $param -fmt 6 1'
'set lev 850'
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
'set lev 850'
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

export param='GP850'
if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson
   export index=$($WGRIB $GRBFILE | grep "HGT"| grep "850 mb" | awk -F ":" '{print $1}')
   export CINT=20
   time $ISOBANDS $GRBFILE -b $index -a $ATTRNAME  -f 'GeoJSON' -i $CINT $OUTPUTPATH/$param.geojson
fi

if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

export param='DIR850'

if [ ! -f $OUTPUTPATH/$param.geojson -o 1$(stat -c %s $OUTPUTPATH/$param.geojson) -le 111111 ] ; then
   [ ! -f $OUTPUTPATH/$param.geojson ] || /bin/rm $OUTPUTPATH/$param.geojson

   export index_u=$($WGRIB $GRBFILE | grep "UGRD"| grep "850 mb" | awk -F ":" '{print $1}')
   export index_v=$($WGRIB $GRBFILE | grep "VGRD"| grep "850 mb" | awk -F ":" '{print $1}')

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
'd skip(1000*ugrd850mb.1+ugrd850mb.2,$WIND_THIN_FACTOR)'
'quit'
EOF

time $GRADS -lbc '$param.gs'
time $OGR2OGR -f 'GeoJSON' $OUTPUTPATH/$param.geojson $param.shp ; sed -i s"/GRID_VALUE/DN/"g $OUTPUTPATH/$param.geojson

fi
if [ ! -f $OUTPUTPATH/$param.mbtiles -o 1$(stat -c %s $OUTPUTPATH/$param.mbtiles) -le 124576 ] ; then
   time $TPCANOE  --layer='MPS_'$param --name='MPS_'$param  -o $OUTPUTPATH/$param.mbtiles $OUTPUTPATH/$param.geojson
fi

