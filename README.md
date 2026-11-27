# ECE 759 Final Project: LDPC Decoding (OpenMP)


OpenMP version of the serial LDPC decoder (`armaan-serial` branch). It decodes many codewords in
parallel and produces bit-identical output to the serial decoder.

## Algorithm

Same decoder as the serial branch:
- (3,6)-regular PEG codes, rate 1/2
- BPSK over AWGN, LLR = 2y/σ²
- normalized min-sum, α = 0.75, at most 20 iterations, early stop when H·ĉᵀ = 0

## Parallelization

Codewords are independent, so the loop over codewords is parallelized (`decode_batch_omp` in `decoder.cpp`):

```cpp
#pragma omp parallel
{
    DecoderWorkspace ws(code);              // private v2c/c2v message arrays per thread
    #pragma omp for schedule(dynamic, 4)    // or schedule(static) with --schedule static
    for (int f = 0; f < frames; ++f)
        iters[f] = decode_frame(code, llr + f*n, bits + f*n, p, ws, &ok);
}
```

- **Private workspace:** each thread has its own message arrays, so there are no data races.
  Each codeword writes only its own output slice.
- **Load balance:** codewords need 0–20 iterations depending on noise. `schedule(dynamic, 4)` lets idle
  threads take more work; `schedule(static)` is included for comparison.
- **Correctness:** each codeword runs the same `decode_frame` as serial, so the output must be identical.
  `./ldpc compare` checks every bit and every iteration count.

## Files

| File | Contents |
|---|---|
| `ldpc.hpp` | Shared data structures. `DecodeParams` adds `threads` and `dynamic_sched` |
| `ldpc_code.cpp` | H matrix construction (PEG), alist I/O, encoder (unchanged from serial) |
| `channel.cpp` | AWGN channel, test data, file I/O (unchanged from serial) |
| `decoder.cpp` | `decode_frame`, `decode_batch_serial`, `decode_batch_omp` |
| `main.cpp` | Command-line program. Adds `--threads` and `--schedule`, and a `threads` CSV column |
| `Makefile` | Builds `./ldpc` with `-fopenmp` |
| `run_euler.sh` | Serial experiments (same as the serial branch) |
| `run_omp_euler.sh` | OpenMP strong-scaling experiment |
| `plot.py` | Plots; adds `results/omp_scaling.png` |

## Build and test

```bash
make                 # Linux / Euler (g++ -fopenmp)
./ldpc test          # must end with "ALL TESTS PASSED"
```

On macOS, Apple clang has no OpenMP. Either build without it (`--threads` then runs serially):
```bash
make OMPFLAGS=
```
or install libomp:
```bash
brew install libomp
make OMPFLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp"
```

## Run on Euler

```bash
sbatch run_omp_euler.sh      # requests 16 cores; output in ldpc_omp-<jobid>.out
squeue -u $USER
```

The script:
1. builds and tests the code
2. generates 4000 codewords (n = 2048, Eb/N0 = 1.75 dB) and decodes them serially (baseline + reference output)
3. decodes the same file with 1, 2, 4, 8, 16 threads, using both `dynamic` and `static` schedules
4. runs `./ldpc compare` against the serial output after every run (each must print `IDENTICAL`)
5. writes `results/omp_scaling_dynamic.csv`, `results/omp_scaling_static.csv` and `results/omp_scaling.png`

## Individual commands

```bash
./ldpc genframes --code codes/peg_2048_3_6.alist --ebno 1.75 --frames 4000 --out data/f.bin
./ldpc decode --code codes/peg_2048_3_6.alist --in data/f.bin --out results/serial.dec
./ldpc decode --code codes/peg_2048_3_6.alist --in data/f.bin --threads 8 --schedule dynamic --out results/omp8.dec
./ldpc compare results/serial.dec results/omp8.dec
./ldpc sweep  --code codes/peg_2048_3_6.alist --sweep 1.0:0.5:2.5 --frames 2000 --threads 8
```

| Option | Meaning | Default |
|---|---|---|
| `--threads N` | number of OpenMP threads (omit for serial) | serial |
| `--schedule` | `dynamic` (chunks of 4) or `static` | dynamic |

The other options (`--iters`, `--alpha`, `--no-early-stop`, `--repeat`, `--seed`) and the output
columns are the same as in the serial branch.

## Results (Euler)

<!-- Fill in from results/omp_scaling_*.csv after running run_omp_euler.sh. Speedup = serial time_ms / time_ms. Efficiency = speedup / threads. -->

Hardware: <!-- CPU model (lscpu), cores, g++ version -->

n = 2048, 4000 codewords, Eb/N0 = 1.75 dB:

| Threads | Time dynamic (ms) | Speedup dynamic | Time static (ms) | Speedup static | Output vs serial |
|---|---|---|---|---|---|
| serial | | 1.00 | | 1.00 | — |
| 1 | | | | | IDENTICAL |
| 2 | | | | | IDENTICAL |
| 4 | | | | | IDENTICAL |
| 8 | | | | | IDENTICAL |
| 16 | | | | | IDENTICAL |

Plot: `results/omp_scaling.png` shows speedup vs thread count for both schedules against ideal linear speedup.

## Notes

- No `-ffast-math`. OpenMP results are bit-identical to serial because each codeword is decoded by
  one thread using the same code path.
- Written with AI assistance (Claude), following the team's Tanner-graph notes.