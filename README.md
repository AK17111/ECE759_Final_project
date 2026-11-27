# ECE 759 Final Project: LDPC Decoding (Serial Baseline)


Serial C++ implementation of an LDPC decoder (normalized min-sum, flooding schedule). It is the
reference that the OpenMP (`parallel-openmp` branch) and CUDA versions are checked against bit for bit.

## Algorithm

- **Code:** (3,6)-regular LDPC codes, rate 1/2, built with Progressive Edge Growth (PEG), no 4-cycles.
  Lengths n = 96, 2048, 4096 and 8192 bits.
- **Encoder:** systematic encoder from Gaussian elimination of H over GF(2).
- **Channel:** BPSK (0 → +1, 1 → −1) over AWGN, with σ² = 1 / (2 · R · Eb/N0). Channel LLR = 2y/σ².
- **Decoder:** normalized min-sum with α = 0.75 and at most 20 iterations. It stops early when every
  parity check is satisfied (H·ĉᵀ = 0).
  - Check update: `c2v = α · Π sign(other v2c) · min |other v2c|`. Uses the two-minimum trick, one pass per check.
  - Bit update: `total = L_ch + Σ c2v`, hard decision `total < 0`, `v2c = total − c2v`.

## Files

| File | Contents |
|---|---|
| `ldpc.hpp` | Shared data structures: `LdpcCode` (H as CSR edge lists), `FrameSet`, `DecodeParams`, `DecoderWorkspace` |
| `ldpc_code.cpp` | H matrix: `code_from_rows`, alist I/O, `make_peg_code`, `count_4cycles`, `build_encoder`, `encode`, `syndrome_ok` |
| `channel.cpp` | AWGN channel and test data: `noise_sigma`, `make_frames`, binary frame and output file I/O |
| `decoder.cpp` | Decoder: `decode_frame` (one codeword), `decode_batch_serial` (all codewords) |
| `main.cpp` | Command-line program (`test`, `gencode`, `genframes`, `sweep`, `decode`, `compare`), timing, BER/FER, unit tests |
| `Makefile` | Builds `./ldpc` |
| `run_euler.sh` | Slurm job that runs every serial experiment and the plots |
| `plot.py` | Turns `results/*.csv` into `results/*.png` |
| `codes/` | Generated H matrices (alist format) |
| `results/` | CSV results and plots |

## Build and test

```bash
make
./ldpc test        # 12 checks, must end with "ALL TESTS PASSED"
```

## Run on Euler

```bash
sbatch run_euler.sh          # ~3 min; output in ldpc_serial-<jobid>.out, results/ and plots
squeue -u $USER
```

The script builds and tests the code, generates the four codes, and runs:
1. BER/FER/throughput vs Eb/N0 for all code lengths → `results/serial_sweep.csv`
2. early termination on vs off (n = 2048) → `results/early_stop_on.csv`, `results/early_stop_off.csv`
3. α = 1.0 vs α = 0.75 (n = 2048) → `results/alpha_1.00.csv`, `results/alpha_0.75.csv`
4. the shared test file `data/n2048_1.75dB.bin` and its serial output, for cross-checking the parallel versions

## Individual commands

```bash
./ldpc gencode   --n 2048 --out codes/peg_2048_3_6.alist
./ldpc sweep     --code codes/peg_2048_3_6.alist --sweep 1.0:0.25:2.5 --frames 2000 --csv results/x.csv
./ldpc genframes --code codes/peg_2048_3_6.alist --ebno 1.75 --frames 1000 --out data/f.bin
./ldpc decode    --code codes/peg_2048_3_6.alist --in data/f.bin --out results/serial.dec
./ldpc compare   results/serial.dec results/other.dec
```

| Option | Meaning | Default |
|---|---|---|
| `--iters` | maximum iterations | 20 |
| `--alpha` | normalization factor | 0.75 |
| `--no-early-stop` | always run every iteration | off |
| `--repeat` | time N runs, keep the fastest | 1 |
| `--seed` | random seed for the test data | 7 |

## Output columns

| Column | Meaning |
|---|---|
| `EbN0_dB` | channel noise level (higher = less noise) |
| `frames` | number of codewords decoded |
| `BER` | wrong payload bits / total payload bits |
| `FER` | codewords with ≥1 payload error / total codewords |
| `avg_iter` | mean decoder iterations per codeword |
| `time_ms` | wall time of decoding only (test-data generation excluded) |
| `info_Mbps`, `coded_Mbps` | decoded payload / coded bits per second (10⁶ bits/s) |

## Tests (`./ldpc test`)

- Hand-worked 6-bit example: one misread bit is corrected in one iteration, with the expected check messages.
- Single check node matches the worked example in the team's Tanner-graph notes (+3.2, −1.1, −4.5 → +1.1).
- The PEG code is (3,6)-regular with zero 4-cycles.
- Encoder output satisfies H·cᵀ = 0 for 200 random payloads.
- Noise-free input decodes in 0 iterations; 10 injected low-confidence errors per codeword are all corrected.
- alist and frame-file round trips.

## Results (Euler)

<!-- Fill in after running run_euler.sh on Euler: CPU model (lscpu), g++ version, and numbers from results/*.csv -->

Hardware: <!-- e.g. Intel Xeon ..., g++ x.y, -O3, 1 thread -->

| n | Eb/N0 (dB) | FER | BER | avg iterations | payload Mb/s |
|---|---|---|---|---|---|
| 2048 | 1.5 | | | | |
| 2048 | 2.0 | | | | |
| 4096 | 2.0 | | | | |
| 8192 | 2.0 | | | | |

Plots:
- `results/ber_fer.png`: BER and FER vs Eb/N0 for each code length
- `results/iters_throughput.png`: iterations and throughput vs Eb/N0
- `results/early_stop.png`: decode time with and without early termination
- `results/alpha.png`: FER for α = 1.0 vs 0.75

## Notes

- No `-ffast-math`, and a fixed summation order (`var_edge`), so floating-point results are reproducible
  and the parallel versions can match the serial output exactly.
- Test data uses fixed seeds. Numbers can differ slightly between macOS and Linux because their standard
  libraries generate normal random numbers differently.
- Written with AI assistance (Claude), following the team's Tanner-graph notes.