"""plot.py - turn results/*.csv into results/*.png (needs matplotlib)."""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def load(name):
    with open(os.path.join(RES, name)) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def by_n(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[int(r["n"])].append(r)
    return dict(sorted(groups.items()))


sweep = by_n(load("serial_sweep.csv"))

# 1) BER and FER vs Eb/N0
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for n, rows in sweep.items():
    xs = [r["ebno_db"] for r in rows]
    for ax, key in zip(axes, ["ber", "fer"]):
        pts = [(x, r[key]) for x, r in zip(xs, rows) if r[key] > 0]
        ax.semilogy([p[0] for p in pts], [p[1] for p in pts], "o-", label=f"n = {n}")
for ax, title in zip(axes, ["Bit error rate (payload bits)", "Frame error rate"]):
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
fig.suptitle("Serial normalized min-sum (alpha = 0.75, max 20 iterations), (3,6)-regular PEG codes")
fig.tight_layout()
fig.savefig(os.path.join(RES, "ber_fer.png"), dpi=130)

# 2) Average iterations and throughput vs Eb/N0
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for n, rows in sweep.items():
    xs = [r["ebno_db"] for r in rows]
    axes[0].plot(xs, [r["avg_iters"] for r in rows], "o-", label=f"n = {n}")
    axes[1].plot(xs, [r["info_mbps"] for r in rows], "o-", label=f"n = {n}")
axes[0].set_title("Average iterations (early stopping on)")
axes[1].set_title("Serial throughput (payload Mb/s)")
for ax in axes:
    ax.set_xlabel("Eb/N0 (dB)")
    ax.grid(True, alpha=0.3)
    ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RES, "iters_throughput.png"), dpi=130)

# 3) Early stop on vs off
on, off = load("early_stop_on.csv"), load("early_stop_off.csv")
fig, ax = plt.subplots(figsize=(6, 4.2))
ax.plot([r["ebno_db"] for r in on], [r["time_ms"] for r in on], "o-", label="early stop on")
ax.plot([r["ebno_db"] for r in off], [r["time_ms"] for r in off], "s--", label="fixed 20 iterations")
ax.set_xlabel("Eb/N0 (dB)")
ax.set_ylabel("Decode time for 1000 frames (ms)")
ax.set_title("Early termination, n = 2048")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RES, "early_stop.png"), dpi=130)

# 4) alpha comparison
fig, ax = plt.subplots(figsize=(6, 4.2))
for name, label in [("alpha_1.00.csv", "plain min-sum (alpha = 1.0)"),
                    ("alpha_0.75.csv", "normalized (alpha = 0.75)")]:
    rows = [r for r in load(name) if r["fer"] > 0]
    ax.semilogy([r["ebno_db"] for r in rows], [r["fer"] for r in rows], "o-", label=label)
ax.set_xlabel("Eb/N0 (dB)")
ax.set_ylabel("Frame error rate")
ax.set_title("Effect of the 0.75 normalization, n = 2048")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RES, "alpha.png"), dpi=130)
print("plots written to results/")