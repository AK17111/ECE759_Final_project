#!/usr/bin/env bash
# OpenMP scaling experiment. On Euler:  sbatch run_omp_euler.sh     Anywhere:  bash run_omp_euler.sh
# Check the partition name (-p) against your course's Euler instructions.
#SBATCH -J ldpc_omp
#SBATCH -p instruction
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 0-00:30:00
#SBATCH -o ldpc_omp-%j.out
#SBATCH -e ldpc_omp-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load gcc 2>/dev/null || true

CORES=${SLURM_CPUS_PER_TASK:-$(nproc)}
echo "cores available: $CORES"
mkdir -p codes results data
make
./ldpc test

CODE=codes/peg_2048_3_6.alist
[ -f $CODE ] || ./ldpc gencode --n 2048 --seed 1 --out $CODE
FRAMES=data/n2048_1.75dB_4000.bin
[ -f $FRAMES ] || ./ldpc genframes --code $CODE --ebno 1.75 --frames 4000 --out $FRAMES

# Serial baseline (threads column = 0) + golden output
rm -f results/omp_scaling.csv
./ldpc decode --code $CODE --in $FRAMES --repeat 3 --out results/serial_4000.dec --csv results/omp_scaling.csv

# 1) Strong scaling: same 4000 codewords, more and more threads. Dynamic vs static schedule.
for sched in dynamic static; do
    t=1
    while [ $t -le $CORES ]; do
        ./ldpc decode --code $CODE --in $FRAMES --repeat 3 --threads $t --schedule $sched \
            --out results/omp_${t}_${sched}.dec --csv results/omp_scaling_${sched}.csv.tmp
        ./ldpc compare results/serial_4000.dec results/omp_${t}_${sched}.dec   # must say IDENTICAL
        t=$((t * 2))
    done
done
# combine (serial row first, then dynamic, then static)
{ cat results/omp_scaling.csv; tail -n +2 results/omp_scaling_dynamic.csv.tmp; } > results/omp_scaling_dynamic.csv
{ cat results/omp_scaling.csv; tail -n +2 results/omp_scaling_static.csv.tmp;  } > results/omp_scaling_static.csv
rm -f results/*.csv.tmp results/omp_*_*.dec

python3 plot.py || echo "plotting skipped (install matplotlib: pip3 install --user matplotlib)"