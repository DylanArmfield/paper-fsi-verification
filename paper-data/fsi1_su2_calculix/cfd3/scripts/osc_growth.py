#!/usr/bin/env python3
"""Growth or decay of the lift oscillation in an SU2 CFD3 run (history*.csv, CL = lift in N/m).

High-pass the lift (minus a moving average over ~one shedding period), then over consecutive
windows of --win seconds report the oscillation half-range and its dominant frequency. The growth
rate sigma (1/s) is a least-squares fit of log(half-range) against time over windows after --skip s
from the start of the file set; sigma < 0 means the shedding mode is damped.

Usage: osc_growth.py CASEDIR|HISTORY.csv [...] [--win 0.5] [--skip 0.5]
"""
import argparse
import csv
import glob
import re

import numpy as np


def read(case):
    """t, lift, drag, inner iterations from a SU2 case folder or a published histories/<run>.csv."""
    if case.endswith(".csv"):
        R = list(csv.DictReader(open(case)))
        return tuple(np.array([float(r[k]) for r in R]) for k in ("t", "lift", "drag", "inner_iters"))
    d = {}
    for fn in sorted(glob.glob(f"{case}/history*.csv")):
        for r in csv.DictReader(open(fn)):
            r = {k.strip().strip('"'): v.strip() for k, v in r.items()}
            d[int(r["Time_Iter"])] = (float(r["CL"]), float(r["CD"]), int(r["Inner_Iter"]) + 1)
    n = sorted(d)
    dt = float(re.search(r"(?m)^TIME_STEP=\s*(\S+)", open(f"{case}/cfd3.cfg").read()).group(1))
    return (np.array([(k + 1) * dt for k in n]), np.array([d[k][0] for k in n]), np.array([d[k][1] for k in n]),
            np.array([d[k][2] for k in n]))


def growth(t, L, it, win=0.5, skip=0.5):
    """Half-range and dominant frequency of the high-passed lift per window; growth rate sigma (1/s)."""
    dt = t[1] - t[0]
    w = max(3, int(round(0.227 / dt)) | 1)
    hp = L - np.convolve(L, np.ones(w) / w, mode="same")
    rows = []
    x = t[0] + skip
    while x + win <= t[-1] - 0.5 * w * dt:
        s = (t >= x) & (t < x + win)
        seg = hp[s]
        f = np.fft.rfftfreq(len(seg), dt)
        A = np.abs(np.fft.rfft(seg - seg.mean()))
        rows.append((x, 0.5 * (seg.max() - seg.min()), f[1 + A[1:].argmax()], it[s].mean(), it[s].max()))
        x += win
    tt = np.array([r[0] + 0.5 * win for r in rows]); amp = np.array([r[1] for r in rows])
    sig = np.polyfit(tt, np.log(amp), 1)[0] if len(rows) > 1 else float("nan")
    return sig, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--win", type=float, default=0.5)
    ap.add_argument("--skip", type=float, default=0.5)
    a = ap.parse_args()
    for c in a.cases:
        t, L, D, it = read(c)
        sig, rows = growth(t, L, it, a.win, a.skip)
        print(f"=== {c}: t {t[0]:g}-{t[-1]:g} s, lift {L[-1]:.4f}, drag {D[-1]:.4f} at end; "
              f"growth rate sigma {sig:+.3f} 1/s (x{np.exp(sig * 0.2275):.3f} per shedding period)")
        for r in rows:
            print(f"    t {r[0]:6.2f}-{r[0] + a.win:6.2f}: half-range {r[1]:.3e}  f {r[2]:5.2f} Hz  inner its mean {r[3]:.1f} max {r[4]}")

if __name__ == "__main__":
    main()
