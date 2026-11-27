"""plot.py - turn results/*.csv into results/*.png (needs matplotlib). Each plot is made only if its CSV exists."""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def have(name):
    return os.path.exists(os.path.join(RES, name))


def load(name):
    with open(os.path.join(RES, name)) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(RES, name), dpi=130)
    print("wrote results/" + name)


# 1) BER and FER vs Eb/N0, and 2) iterations + throughput vs Eb/N0
if have("serial_sweep.csv"):
    sweep = defaultdict(list)
    for r in load("serial_sweep.csv"):
        sweep[int(r["n"])].append(r)
    sweep = dict(sorted(sweep.items()))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for n, rows in sweep.items():
        for ax, key in zip(axes, ["ber", "fer"]):
            pts = [(r["ebno_db"], r[key]) for r in rows if r[key] > 0]   # log scale cannot show 0
            ax.semilogy([p[0] for p in pts], [p[1] for p in pts], "o-", label=f"n = {n}")
    for ax, title in zip(axes, ["Bit error rate (payload bits)", "Frame error rate (codewords)"]):
        ax.set_xlabel("Eb/N0 (dB)  - higher = less noise")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    fig.suptitle("Serial normalized min-sum (alpha = 0.75, max 20 iterations), (3,6)-regular PEG codes")
    save(fig, "ber_fer.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for n, rows in sweep.items():
        xs = [r["ebno_db"] for r in rows]
        axes[0].plot(xs, [r["avg_iters"] for r in rows], "o-", label=f"n = {n}")
        axes[1].plot(xs, [r["info_mbps"] for r in rows], "o-", label=f"n = {n}")
    axes[0].set_title("Average iterations (early stopping on)")
    axes[0].set_ylabel("iterations")
    axes[1].set_title("Serial throughput")
    axes[1].set_ylabel("payload Mb/s")
    for ax in axes:
        ax.set_xlabel("Eb/N0 (dB)")
        ax.grid(True, alpha=0.3)
        ax.legend()
    save(fig, "iters_throughput.png")

# 3) Early stop on vs off
if have("early_stop_on.csv") and have("early_stop_off.csv"):
    on, off = load("early_stop_on.csv"), load("early_stop_off.csv")
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot([r["ebno_db"] for r in on], [r["time_ms"] for r in on], "o-", label="early stop on")
    ax.plot([r["ebno_db"] for r in off], [r["time_ms"] for r in off], "s--", label="fixed 20 iterations")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("Decode time for 1000 codewords (ms)")
    ax.set_title("Early termination, n = 2048")
    ax.grid(True, alpha=0.3)
    ax.legend()
    save(fig, "early_stop.png")

# 4) alpha comparison
if have("alpha_1.00.csv") and have("alpha_0.75.csv"):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for name, label in [("alpha_1.00.csv", "plain min-sum (alpha = 1.0)"),
                        ("alpha_0.75.csv", "scaled (alpha = 0.75)")]:
        rows = [r for r in load(name) if r["fer"] > 0]
        ax.semilogy([r["ebno_db"] for r in rows], [r["fer"] for r in rows], "o-", label=label)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("Frame error rate")
    ax.set_title("Effect of the 0.75 scaling, n = 2048")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save(fig, "alpha.png")

# 5) OpenMP strong scaling
if have("omp_scaling_dynamic.csv"):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    max_t = 1
    for sched, style in [("dynamic", "o-"), ("static", "s--")]:
        if not have(f"omp_scaling_{sched}.csv"):
            continue
        rows = load(f"omp_scaling_{sched}.csv")
        serial_ms = [r for r in rows if r["threads"] == 0][0]["time_ms"]
        par = [r for r in rows if r["threads"] > 0]
        ts = [r["threads"] for r in par]
        max_t = max(max_t, max(ts))
        ax.plot(ts, [serial_ms / r["time_ms"] for r in par], style, label=f"schedule({sched})")
    ax.plot([1, max_t], [1, max_t], ":", color="gray", label="ideal (speedup = threads)")
    ax.set_xlabel("OpenMP threads")
    ax.set_ylabel("Speedup over serial (x)")
    ax.set_title("OpenMP strong scaling, n = 2048, 4000 codewords, 1.75 dB")
    ax.grid(True, alpha=0.3)
    ax.legend()
    save(fig, "omp_scaling.png")