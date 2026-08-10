#!/bin/bash
#SBATCH --error=/scratch/lus/arw/model/log/mpas2wrf.%a.out
#SBATCH --output=/scratch/lus/arw/model/log/mpas2wrf.%a.out
#SBATCH --job-name=MPAS2WRF
#SBATCH --requeue
#SBATCH --partition=opr
#SBATCH --exclusive
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=0-228:1
#SBATCH --export=ALL


# -----------------------------
# Environment
# -----------------------------
. /scratch/lus/arw/model/slurm/modelenv.sh

set -x

# Keep stdout/stderr for each array index
exec >>/scratch/lus/arw/model/log/mpas2wrf.${SLURM_ARRAY_TASK_ID}.out 2>&1

ulimit -s unlimited
ulimit -c 0

echo "============================================================"
echo "Job      : ${SLURM_JOB_ID:-NA}"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Host     : $(hostname)"
echo "Date     : $(date -Is)"
echo "============================================================"

# metgrid is MPI; avoid accidental threading overhead
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# -----------------------------
# Compiler / modules
# -----------------------------
export COMPILER=intel

module purge
unset LD_LIBRARY_PATH
unset LIBRARY_PATH
unset CPATH
unset CMAKE_PREFIX_PATH
unset PKG_CONFIG_PATH

case "$COMPILER" in
  cce)   module load PrgEnv-cray/8.5.0 ;;
  intel) module load PrgEnv-intel/8.5.0 ;;
  aocc)  module load PrgEnv-aocc/8.5.0 ;;
  nvhpc) module load PrgEnv-nvhpc/8.5.0 ;;
  *)     echo "Unknown COMPILER=$COMPILER"; exit 2 ;;
esac

# I/O libs (Cray)
module load cray-netcdf/4.9.0.13
module load cray-parallel-netcdf/1.12.3.13
module load craype-x86-rome

if [[ "${USE_UCX:-no}" == "yes" ]]; then
  module swap craype-network-ofi craype-network-ucx
  module swap cray-mpich cray-mpich-ucx
fi

# -----------------------------
# Paths
# -----------------------------
export MPAS2WRF="$MODEL_ROOT/src/util/convert.$COMPILER/metgrid/src/metgrid.exe"
export EXE=${MPAS2WRF}

# Work in fast local scratch
cd "${TMPDIR:?TMPDIR not set}"

# -----------------------------
# Dates
# -----------------------------
export CYCLE="$(echo "$START_DATE" | cut -c9-10)"

rm -f "$MODEL_RUN/r$CYCLE/OK.MPAS2WRF.$START_DATE" || true

export START_YEAR="$(echo "$START_DATE" | cut -c1-4)"
export START_MONTH="$(echo "$START_DATE" | cut -c5-6)"
export START_DAY="$(echo "$START_DATE" | cut -c7-8)"
export START_HOUR="$(echo "$START_DATE" | cut -c9-10)"

export ACTUAL_DATE="$($SMSDATE +$SLURM_ARRAY_TASK_ID $START_DATE)"
export ACTUAL_YEAR="$(echo "$ACTUAL_DATE" | cut -c1-4)"
export ACTUAL_MONTH="$(echo "$ACTUAL_DATE" | cut -c5-6)"
export ACTUAL_DAY="$(echo "$ACTUAL_DATE" | cut -c7-8)"
export ACTUAL_HOUR="$(echo "$ACTUAL_DATE" | cut -c9-10)"

# -----------------------------
# MPAS grid selection
# -----------------------------
case "${MPAS_RESOLUTION:?MPAS_RESOLUTION not set}" in
  10)    XNN=x1;  NUM=5898242; DT=60.0; DX=10000.0; UPWIND=0.0; DUST="${DUST:-}";;
  12)    XNN=x1;  NUM=4096002; DT=72.0; DX=12000.0; UPWIND=0.0; DUST="${DUST:-}";;
  15)    XNN=x1;  NUM=2621442; DT=90.0; DX=15000.0; UPWIND=0.0; DUST="${DUST:-}";;
  15_3)  XNN=x5;  NUM=8060930; DT=20.0; DX=3000.0 ; UPWIND=0.0; DUST="${DUST:-}";;
  15_5)  XNN=x3;  NUM=4763609; DT=36.0; DX=5000.0 ; UPWIND=0.0; DUST="${DUST:-}";;
  30_5)  XNN=x6;  NUM=2819097; DT=36.0; DX=5000.0 ; UPWIND=0.0; DUST="${DUST:-}";;
  60_10) XNN=x6;  NUM=999426;  DT=60.0; DX=10000.0; UPWIND=0.5; DUST="${DUST:-}";;
  60_3)  XNN=x20; NUM=835586;  DT=20.0; DX=3000.0 ; UPWIND=0.5; DUST="${DUST:-}";;
  46_12) XNN=x4;  NUM=655362;  DT=72.0; DX=12000.0; UPWIND=0.5; DUST="${DUST:-}";;
  60_15) XNN=x4;  NUM=535554;  DT=90.0; DX=15000.0; UPWIND=0.5; DUST="${DUST:-}";;
  *) echo "Unsupported MPAS_RESOLUTION=$MPAS_RESOLUTION"; exit 0;;
esac
export MPAS_GRID="$XNN.$NUM.static$DUST.nc"

# -----------------------------
# I/O tuning (Cray Shasta + Lustre, 8 OSTs)
# -----------------------------
# PnetCDF: disable debug/safe checks for speed
export PNETCDF_SAFE_MODE=0

# Your filesystem has 8 OSTs -> do NOT exceed 8 stripes.
export STRIPE=8
export STRIPE_UNIT=$((8*1024*1024))   # 8 MiB (good default for large netcdf writes)
export CB_NODES=4                     # try 2 if you want to test further

# Apply striping to output directory (best practice; ignore failure if not Lustre or no perms)
OUTDIR="$MODEL_RUN/r$CYCLE"
if command -v lfs >/dev/null 2>&1; then
  lfs setstripe -c "$STRIPE" -S "${STRIPE_UNIT}" "$OUTDIR" >/dev/null 2>&1 || true
fi

# MPI-IO hints: make sure the pattern matches actual output file path/name.
# metgrid writes "met_em*" files (in your case netcdf). Adjust if your filenames differ.
export MPICH_MPIIO_HINTS="\
$OUTDIR/met_em*:striping_factor=$STRIPE:striping_unit=$STRIPE_UNIT:cb_nodes=$CB_NODES:cray_cb_nodes_multiplier=$CB_NODES:cray_cb_write_lock_mode=2"

# Optional debugging while tuning (comment out for production)
# export MPICH_MPIIO_HINTS_DISPLAY=1
# export MPICH_MPIIO_STATS=1

# Coll tuning: keep defaults unless you measured benefit.
# export MPICH_RANK_REORDER_METHOD=0

# -----------------------------
# Build namelist
# -----------------------------
sed -e "s/BASE_DATE/${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}/g" \
    -e "s/YEAR/${ACTUAL_YEAR}/g" \
    -e "s/MONTH/${ACTUAL_MONTH}/g" \
    -e "s/DAY/${ACTUAL_DAY}/g" \
    -e "s/HOUR/${ACTUAL_HOUR}/g" \
    -e "s/DOMAIN/${MPAS_DOMAIN}/g" \
    -e "s/CYCLE/${CYCLE}/g" \
    -e "s/VERSION/${WRF_VERSION}/g" \
    -e "s/COMPILER/${COMPILER}/g" \
    -e "s/LBC_FREQUENCY/3600/g" \
    -e "s/MPAS_GRID/${MPAS_GRID}/g" \
    -e "s/XNN/${XNN}/g" \
    -e "s/NUM/${NUM}/g" \
    -e "s/DUST/${DUST}/g" \
    -e "s/COMPILER/${COMPILER}/g" \
    "$MODEL_ROOT/name/wrf/$WRF_VERSION/namelist.metgrid.mpas2wrf" > namelist.wps

echo "----- namelist.wps -----"
cat namelist.wps
echo "------------------------"

# Logs & links
ln -sf "$MODEL_LOG/mpas2wrf.$(printf "%03d" "$SLURM_ARRAY_TASK_ID").log" metgrid.log.0000 || true
ln -sf "$MESH_DIR/$MPAS_GRID" . || true

# -----------------------------
# Wait for target MPAS file (avoid du in a loop; use stat only)
# -----------------------------
TARGET="$MODEL_RUN/r$CYCLE/MPAS.${ACTUAL_YEAR}-${ACTUAL_MONTH}-${ACTUAL_DAY}_${ACTUAL_HOUR}.nc"
MIN_BYTES=$((20000000 * 1024))  # ~20GB (you used du -sk > 20000000)

echo "Waiting for: $TARGET"
for ((i=0; i<180; i++)); do
  if [[ -f "$TARGET" ]]; then
    sz=$(stat -c %s "$TARGET" 2>/dev/null || echo 0)
    if (( sz > MIN_BYTES )); then
      ln -sf "$TARGET" .
      echo "File ready ($sz bytes): $TARGET"
      break
    fi
  fi
  sleep 30
done

if (( i == 180 )); then
  echo "File not found or too small after 90 minutes. Exiting."
  exit 0
fi

# If already processed, skip
OUTFILE="$MODEL_RUN/r$CYCLE/mpsout.d01.${ACTUAL_YEAR}-${ACTUAL_MONTH}-${ACTUAL_DAY}_${ACTUAL_HOUR}_00_00.nc"
if [[ -f "$OUTFILE" ]]; then
  outsz=$(stat -c %s "$OUTFILE" 2>/dev/null || echo 0)
  if [[ -n "${MPASOUT_FILE_SIZE:-}" ]]; then
    if (( outsz > MPASOUT_FILE_SIZE )); then
      echo "Output already exists and is large enough ($outsz > $MPASOUT_FILE_SIZE). Exiting."
      exit 0
    fi
  fi
fi

# -----------------------------
# Run
# -----------------------------
echo "Running metgrid: $EXE"
echo "srun -n $SLURM_NTASKS --cpu-bind=cores $EXE"

$STIME srun -v -l -n "$SLURM_NTASKS" -N "$SLURM_NNODES" --cpu-bind=cores "$EXE"

echo "Done: $(date -Is)"

