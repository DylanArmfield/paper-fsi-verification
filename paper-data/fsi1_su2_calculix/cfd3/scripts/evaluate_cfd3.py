#!/usr/bin/env python3
"""Evaluate SU2 CFD3 runs (history.csv, CD/CL = drag/lift in N/m on cylinder + plate).

Turek-Hron convention: over the last N complete lift periods (lift maxima located by parabolic
interpolation), mean = (max+min)/2, amplitude = (max-min)/2, frequency = 1/period of the lift.
Drag and lift max/min also use parabolic interpolation around the extreme sample. Periodicity check:
the same quantities over each single period of the window; the spread (max-min over the periods) is
reported. Also: inner iterations per step (mean, max, steps at the cap) over the window.

Usage: evaluate_cfd3.py CASEDIR|HISTORY.csv [...] [--periods 4] [--cap 1000] [--csv OUT]
"""
import argparse
import csv
import glob
import math
import re

REF = {"drag_mean": 439.45, "drag_amp": 5.6183, "lift_mean": -11.893, "lift_amp": 437.81, "freq": 4.3956}
# solids4foam CFD3, matched dt sequence: ../../../hronturek_cfd3_runs.csv rows at dt 1e-3 / 5e-4 / 2.5e-4
S4F = {1: {"drag_mean": 446.0105, "drag_amp": 4.0616, "lift_mean": -67.9215, "lift_amp": 296.0431},
       2: {"drag_mean": 441.8067, "drag_amp": 5.5267, "lift_mean": -18.9754, "lift_amp": 424.4350},
       4: {"drag_mean": 440.1482, "drag_amp": 5.7317, "lift_mean": -13.8417, "lift_amp": 438.0441}}


def read(case):
    """t, drag, lift, inner iterations from a SU2 case folder (history*.csv; t = (Time_Iter + 1) dt) or
    from a published histories/<run>.csv (t, drag, lift, inner_iters, ...)."""
    if case.endswith(".csv"):
        R = list(csv.DictReader(open(case)))
        return ([float(r["t"]) for r in R], [float(r["drag"]) for r in R], [float(r["lift"]) for r in R],
                [int(r["inner_iters"]) for r in R])
    d = {}
    for fn in sorted(glob.glob(f"{case}/history*.csv")):
        for r in csv.DictReader(open(fn)):
            r = {k.strip().strip('"'): v.strip() for k, v in r.items()}
            d[int(r["Time_Iter"])] = (float(r["CD"]), float(r["CL"]), int(r["Inner_Iter"]) + 1)
    n = sorted(d)
    dt = float(re.search(r"(?m)^TIME_STEP=\s*(\S+)", open(f"{case}/cfd3.cfg").read()).group(1))
    return [(k + 1) * dt for k in n], [d[k][0] for k in n], [d[k][1] for k in n], [d[k][2] for k in n]


def parab(t, y, i, sign):
    """Extremum of the parabola through samples i-1, i, i+1 (sign +1 max, -1 min)."""
    if i <= 0 or i >= len(y) - 1:
        return t[i], y[i]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    den = y0 - 2 * y1 + y2
    if den == 0:
        return t[i], y1
    s = 0.5 * (y0 - y2) / den
    h = t[i + 1] - t[i]
    return t[i] + s * h, y1 - 0.25 * (y0 - y2) * s


def maxima(t, y):
    return [parab(t, y, i, 1) for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1]]


def stats(t, y, a, b):
    idx = [i for i in range(len(t)) if a <= t[i] <= b]
    imax = max(idx, key=lambda i: y[i]); imin = min(idx, key=lambda i: y[i])
    mx = parab(t, y, imax, 1)[1]; mn = parab(t, y, imin, -1)[1]
    return 0.5 * (mx + mn), 0.5 * (mx - mn)


def evaluate(case, N, cap):
    t, D, Lf, it = read(case)
    pk = maxima(t, Lf)
    pk = [p for p in pk if p[1] > 0.5 * max(Lf[len(Lf) // 2:])] if len(Lf) > 10 else pk
    if len(pk) < N + 1:
        return {"case": case, "t_end": t[-1], "error": f"only {len(pk)} lift maxima"}
    w = pk[-(N + 1):]
    a, b = w[0][0], w[-1][0]
    per = [(w[k][0], w[k + 1][0]) for k in range(N)]
    out = {"case": case, "t_end": t[-1], "window": (a, b), "freq": N / (b - a)}
    for name, y in (("drag", D), ("lift", Lf)):
        m, amp = stats(t, y, a, b)
        out[f"{name}_mean"], out[f"{name}_amp"] = m, amp
        single = [stats(t, y, p, q) for p, q in per]
        out[f"{name}_spread_mean"] = max(s[0] for s in single) - min(s[0] for s in single)
        out[f"{name}_spread_amp"] = max(s[1] for s in single) - min(s[1] for s in single)
    fr = [1 / (q - p) for p, q in per]
    out["freq_spread"] = max(fr) - min(fr)
    wi = [n for x, n in zip(t, it) if a <= x <= b]
    out["inner_mean"], out["inner_max"], out["inner_capped"] = sum(wi) / len(wi), max(wi), sum(n >= cap for n in wi)
    out["inner_mean_all"] = sum(it) / len(it)
    out["steps"] = len(t)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--periods", type=int, default=4)
    ap.add_argument("--cap", type=int, default=1000)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    res = [evaluate(c, a.periods, a.cap) for c in a.cases]
    keys = ["drag_mean", "drag_amp", "lift_mean", "lift_amp", "freq"]
    for r in res:
        print(f"=== {r['case']}  (t_end {r['t_end']:g} s)")
        if "error" in r:
            print("   ", r["error"]); continue
        print(f"    window {r['window'][0]:.4f}-{r['window'][1]:.4f} s ({a.periods} lift periods); "
              f"inner its mean {r['inner_mean']:.1f} (whole run {r['inner_mean_all']:.1f}), max {r['inner_max']}, "
              f"at cap {r['inner_capped']}")
        for k in keys:
            sp = r.get(k.replace("_mean", "_spread_mean").replace("_amp", "_spread_amp")) if k != "freq" else r["freq_spread"]
            print(f"    {k:10s} {r[k]:12.5f}   vs Featflow {100*(r[k]-REF[k])/abs(REF[k]):+7.2f} %   "
                  f"period-to-period spread {sp:.2e}")
    if a.csv:
        cols = ["case", "t_end", "steps"] + keys + ["drag_spread_mean", "drag_spread_amp", "lift_spread_mean",
                                                    "lift_spread_amp", "freq_spread", "inner_mean", "inner_max",
                                                    "inner_capped", "inner_mean_all"]
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(cols)
            for r in res:
                if "error" not in r:
                    w.writerow([r[c] for c in cols])


if __name__ == "__main__":
    main()
