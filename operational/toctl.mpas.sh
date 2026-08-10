#!/bin/bash
#SBATCH --error=/scratch/lus/arw/model/log/toctl.mpas.out
#SBATCH --output=/scratch/lus/arw/model/log/toctl.mpas.out
#SBATCH --job-name=CTL.MPS
#SBATCH --partition=web
#SBATCH --nodes=1
#SBATCH --time=1:45:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL

> /scratch/lus/arw/model/log/toctl.mpas.out

cd $TMPDIR

. $HOME/model/slurm/modelenv.sh

set -x

export LOC_GRIB=/scratch/lus/data/LOC
export RANGE=229

export CYCLE=${START_DATE:8:2}
export START_YEAR=${START_DATE:0:4}
export START_MONTH=${START_DATE:4:2}
export START_DAY=${START_DATE:6:2}

export YYDATE=$($SMSDATE -24 $START_DATE)
export YDATE=${YYDATE:0:8}
export XDATE=${START_DATE:0:8}

D="$LOC_GRIB/$START_DATE"

# Generate GRADS control file
$GRIB2CTL -verf $D/mpsprs.$START_YEAR-$START_MONTH-${START_DAY}_${CYCLE}_00_00.grb1 > GRADS.MPS.CTL

awk -v D="$D" 'NR==1 {print "dset " D "/mpsprs.%y4-%m2-%d2_%h2_00_00.grb1"} NR==2 {print "index " D "/mpsprs.idx"} NR==4 {print "title MPS"} NR==5 {print "options template"} NR!=1 && NR!=2 && NR!=4 && NR!=5 {print}' GRADS.MPS.CTL > tempo

N=$(grep -n tdef tempo | cut -d: -f1 | sed 's/^0*//')

index=$((10#$START_MONTH - 1))

TDEF="tdef $RANGE linear ${CYCLE}Z${START_DAY}${MONTHS[$index]}${START_YEAR} 1hr"
awk -v N="$N" -v TDEF="$TDEF" 'NR==N {$0=TDEF} {print}' tempo > $D/GRADS.MPS.CTL

# Generate GRIB index
$GRIBMAP -big -i $D/GRADS.MPS.CTL

# Update last run
echo $START_DATE > $WEB_DIR/js/MPS.date_of_last_run

