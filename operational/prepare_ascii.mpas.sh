#! /bin/bash
# This script creates empty ascii data files
# Base on the actual date and time (Date and time taken from the system clock)

#SBATCH --error=/scratch/lus/arw/model/log/prepare_ascii.mpas.out
#SBATCH --output=/scratch/lus/arw/model/log/prepare_ascii.mpas.out
#SBATCH --job-name=PREP.ASCII --no-requeue
#SBATCH --ntasks=1
#SBATCH --partition=opr  --no-requeue
#SBATCH --time=0:30:00
#SBATCH --propagate=STACK,CORE
#SBATCH --export=ALL

. $HOME/model/slurm/modelenv.sh

export CYCLE=$(expr $START_DATE | cut -c9-10)
cp $MODEL_ROOT/etc/indexes.stations.mpas indexes
$HOME/model/src/ascii/dev/point.exe -multi \
                                    -nearest \
                                    -d MPS \
                                    -cloc /scratch/lus/arw/model/etc/input_stations_global \
                                    -dummy $CYCLE
