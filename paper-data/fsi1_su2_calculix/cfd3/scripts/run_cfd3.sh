#!/bin/bash
# SU2 CFD3 (rigid flag, unsteady) in one case folder made by scripts/setup_cfd3.py.
# Submit from this cfd3/ folder:  sbatch --ntasks=N --export=ALL,CASE=cases/<case> scripts/run_cfd3.sh
#SBATCH --job-name=cfd3-su2
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
echo "Job $SLURM_JOB_ID ($SLURM_JOB_NAME)  host $(hostname)  tasks ${SLURM_NTASKS:-1}  mem ${SLURM_MEM_PER_NODE:-?} MB  start $(date -Is)"
cd "$SLURM_SUBMIT_DIR/$CASE"
echo "case $PWD"; cat CASE_INFO.txt
set +u; source "${FSI_ENV:-$SLURM_SUBMIT_DIR/../scripts/env-su2-precice.sh}"; set -u   # see ../scripts/env-su2-precice.sh
export OMP_NUM_THREADS=1
date +%s > wall_start.txt
set +e
mpirun -np "$SLURM_NTASKS" SU2_CFD cfd3.cfg > su2.log 2>&1
STATUS=$?
set -e
date +%s > wall_end.txt
echo $STATUS > exit_status.txt
python3 - <<'EOF' || true
import csv, glob
f = sorted(glob.glob("history*.csv"))[-1]
R = [{k.strip().strip('"'): v.strip() for k, v in r.items()} for r in csv.DictReader(open(f))]
it = [int(r["Inner_Iter"]) + 1 for r in R]
print(f"time steps {len(R)}, t_end {float(R[-1]['Cur_Time']):g}, inner its mean {sum(it)/len(it):.1f} max {max(it)}")
EOF
echo "Exit $STATUS  end $(date -Is)  wall $(( $(cat wall_end.txt) - $(cat wall_start.txt) )) s"
exit $STATUS
