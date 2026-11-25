#!/usr/bin/env bash
# Runs every serial experiment. Works two ways:
#   On Euler:   sbatch run_euler.sh
#   Anywhere:   bash run_euler.sh
# Check the partition name (-p) against your course's Euler instructions.
#SBATCH -J ldpc_serial
#SBATCH -p instruction
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-00:30:00
#SBATCH -o ldpc_serial-%j.out
#SBATCH -e ldpc_serial-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load gcc 2>/dev/null || true

mkdir -p codes results data
make
./ldpc test

# 1) The H matrices: (3,6)-regular PEG codes, rate 1/2, no 4-cycles.
for n in 96 2048 4096 8192; do
    [ -f codes/peg_${n}_3_6.alist ] || ./ldpc gencode --n $n --seed 1 --out codes/peg_${n}_3_6.alist
done

# 2) Error rate + speed vs noise (~4 million bits per noise level).
rm -f results/serial_sweep.csv
./ldpc sweep --code codes/peg_96_3_6.alist --sweep 0.5:0.5:4.0 --frames 20000 --repeat 3 --csv results/serial_sweep.csv
for n in 2048 4096 8192; do
    ./ldpc sweep --code codes/peg_${n}_3_6.alist --sweep 0.5:0.25:2.5 --frames $((4096000 / n)) \
        --repeat 3 --csv results/serial_sweep.csv
done

# 3) Early stopping ON vs OFF, n=2048.
rm -f results/early_stop_on.csv results/early_stop_off.csv
./ldpc sweep --code codes/peg_2048_3_6.alist --sweep 1.0:0.25:2.5 --frames 1000 --repeat 3 --csv results/early_stop_on.csv
./ldpc sweep --code codes/peg_2048_3_6.alist --sweep 1.0:0.25:2.5 --frames 1000 --repeat 3 --no-early-stop --csv results/early_stop_off.csv

# 4) Plain min-sum (alpha=1.0) vs scaled (alpha=0.75), n=2048.
rm -f results/alpha_1.00.csv results/alpha_0.75.csv
for al in 1.00 0.75; do
    ./ldpc sweep --code codes/peg_2048_3_6.alist --sweep 1.0:0.25:2.5 --frames 2000 --alpha $al --csv results/alpha_${al}.csv
done

# 5) Shared test file + serial output: OpenMP/CUDA must match this exactly.
./ldpc genframes --code codes/peg_2048_3_6.alist --ebno 1.75 --frames 1000 --out data/n2048_1.75dB.bin
./ldpc decode --code codes/peg_2048_3_6.alist --in data/n2048_1.75dB.bin --out results/serial_n2048_1.75dB.dec

python3 plot.py || echo "plotting skipped (install matplotlib: pip3 install --user matplotlib)"