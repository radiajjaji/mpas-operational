#! /bin/bash
# Script to extract ascii data from grib files

#SBATCH --output=/scratch/lus/arw/model/log/ascii.mpas.%a.out
#SBATCH --error=/scratch/lus/arw/model/log/ascii.mpas.%a.out
#SBATCH --job-name=ASCII.MPAS --no-requeue
#SBATCH --nodes=1 --ntasks-per-node=126 --cpus-per-task=1
#SBATCH --partition=workq  --no-requeue
#SBATCH --time=0:20:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL
#SBATCH --array=0-228:1

. $HOME/model/slurm/modelenv.sh

cd $TMPDIR

#cp $MODEL_ROOT/etc/indexes.stations.mpas indexes

export scdate=$START_DATE
export ech=$SLURM_ARRAY_TASK_ID

export acdate=$(/scratch/lus/dev/bin/smsdate +$ech $scdate)
export ayear=$(expr $acdate | cut -c1-4)
export amonth=$(expr $acdate | cut -c5-6)
export aday=$(expr $acdate | cut -c7-8)
export ahour=$(expr $acdate | cut -c9-10)

$STIME $ASCII_ROOT/dev/point.exe -multi -sound -series  -nearest -d MPS \
        -bdt $START_DATE \
        -fdt ${ayear}${amonth}${aday}${ahour} \
	-cloc /scratch/lus/arw/model/etc/input_stations 

