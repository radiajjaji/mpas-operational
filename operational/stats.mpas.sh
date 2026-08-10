#! /bin/bash
# This script calculates many forecast needed statistics based on  ascii data.
# Base on the actual date and time (Date and time taken from the system clock)
#
#SBATCH --error=/scratch/lus/arw/model/log/stats.mpas.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/stats.mpas.%a.out
#SBATCH --job-name=STATS.MPAS --no-requeue
#SBATCH --nodes=1 --ntasks=126 
#SBATCH --partition=uan,hpclm,slm
#SBATCH --time=0:10:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL
#SBATCH --array=0-80:1

. $HOME/model/slurm/modelenv.sh

cd $TMPDIR

/usr/bin/split -d -l 150 /scratch/lus/arw/model/etc/input_stations
ls -rtl

$HOME/model/src/ascii/dev/point.exe -V -multi -stats \
                                    -bdt $START_DATE \
                                    -nearest \
                                    -d MPS \
                                    -cloc $(pwd)/x$(printf "%02d" $SLURM_ARRAY_TASK_ID)
