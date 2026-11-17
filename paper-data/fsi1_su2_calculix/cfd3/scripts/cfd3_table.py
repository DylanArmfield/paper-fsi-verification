#!/usr/bin/env python3
"""Summary table of the SU2 CFD3 runs from the published histories and logs (cfd3_su2.csv).

For every run: settings (from logs/<run>/CASE_INFO.txt), cores and wall clock (logs/<run>/slurm.out,
walltime.txt), the lift-oscillation growth per shedding period (osc_growth.growth over the whole
history; left empty when the lift oscillation never exceeds 1e-3 N/m, i.e. a steady flow whose
residual fluctuation is solver noise), the first and last 0.5 s half-ranges of the high-passed lift
(lift_osc_first/last), the drag and lift at the end, inner-iteration statistics, and, for the runs that shed
periodically, the Turek-Hron quantities over the last 4 lift periods (evaluate_cfd3.evaluate) with
their deviations from Featflow and from the solids4foam CFD3 run on the same mesh.

Usage (from cfd3/):  python3 scripts/cfd3_table.py [--csv cfd3_su2.csv]
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import evaluate_cfd3 as ev  # noqa: E402
import osc_growth as og  # noqa: E402

RUNS = ["dev-L1-cap40", "diag-L1-b4.1", "diag-L1-b20", "diag-L1-b100", "diag-L1-jst",
        "dev-L2-cap40", "conv-L2-b4.1", "conv-L2-b20"]
SHEDDING = {"dev-L2-cap40", "conv-L2-b4.1", "conv-L2-b20"}
PERIOD = 0.2275  # s, the Featflow shedding period (1 / 4.3956 Hz)


def info(run):
    s = (ROOT / "logs" / run / "CASE_INFO.txt").read_text()
    h = (ROOT / "logs" / run / "slurm.out").read_text().splitlines()[0]
    g = lambda p, f=str: f(re.search(p, s).group(1))
    return dict(run=run, level=g(r"CFD3 level (\d+)", int), cells=g(r"(\d+) cells", int), dt=g(r"dt (\S+) s", float),
                beta=g(r"BETA_FACTOR (\S+?),", float), scheme="JST" if "JST" in s else "FDS",
                inner_cap=g(r"inner iter cap (\d+)", int), inner_minval=g(r"inner minval (\S+?)[;\s]", float),
                restart_from=(re.search(r"restart from (\S+) at", s) or [None, ""])[1],
                job=int(re.search(r"Job (\d+)", h).group(1)), cores=int(re.search(r"tasks (\d+)", h).group(1)),
                wall_s=int((ROOT / "logs" / run / "walltime.txt").read_text()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(ROOT / "cfd3_su2.csv"))
    a = ap.parse_args()
    rows = []
    for run in RUNS:
        r = info(run)
        hist = str(ROOT / "histories" / f"{run}.csv")
        t, L, D, it = og.read(hist)
        sig, win = og.growth(t, L, it)
        hr = [w[1] for w in win]
        steps = len(t)
        r.update(t_start=round(t[0], 4), t_end=round(t[-1], 4), steps=steps, s_per_step=r["wall_s"] / steps,
                 core_hours=r["wall_s"] * r["cores"] / 3600, inner_mean=float(it.mean()), inner_max=int(it.max()),
                 inner_at_cap=int((it >= r["inner_cap"]).sum()),
                 lift_osc_first=hr[0], lift_osc_last=hr[-1],
                 # growth only where there is a shedding oscillation to measure (half-range > 1e-3 N/m)
                 growth_per_period=math.exp(sig * PERIOD) if max(hr) > 1e-3 else float("nan"),
                 drag_end=float(D[-1]), lift_end=float(L[-1]))
        if run in SHEDDING:
            e = ev.evaluate(hist, 4, r["inner_cap"])
            r["window"] = f"{e['window'][0]:.4f}-{e['window'][1]:.4f}"
            s4f = ev.S4F.get(r["level"], {})
            for k in ("drag_mean", "drag_amp", "lift_mean", "lift_amp", "freq"):
                r[k] = e[k]
                r[f"{k}_vs_featflow_pct"] = 100 * (e[k] - ev.REF[k]) / abs(ev.REF[k])
                if k in s4f:
                    r[f"{k}_vs_s4f_same_mesh_pct"] = 100 * (e[k] - s4f[k]) / abs(s4f[k])
            for k in ("drag_spread_mean", "drag_spread_amp", "lift_spread_mean", "lift_spread_amp", "freq_spread"):
                r[k] = e[k]
        rows.append(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(a.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([f"{r[k]:.7g}" if isinstance(r.get(k), float) else r.get(k, "") for k in keys])
    for r in rows:
        s = (f"{r['run']:14s} L{r['level']} beta {r['beta']:5g} {r['scheme']} cap {r['inner_cap']:4d}: "
             f"t {r['t_start']:g}-{r['t_end']:g}, osc {r['lift_osc_first']:.1e}->{r['lift_osc_last']:.1e} "
             f"(x{r['growth_per_period']:.3f}/period), end drag {r['drag_end']:.2f} "
             f"lift {r['lift_end']:.2f}; inner {r['inner_mean']:.1f} ({r['inner_at_cap']} at cap); "
             f"{r['cores']} cores {r['wall_s']} s ({r['s_per_step']:.2f} s/step)")
        if "drag_mean" in r:
            s += (f"\n{'':16s}drag {r['drag_mean']:.2f} +- {r['drag_amp']:.2f}, lift {r['lift_mean']:.2f} +- "
                  f"{r['lift_amp']:.1f}, f {r['freq']:.4f} Hz")
        print(s)
    print(f"wrote {a.csv}")


if __name__ == "__main__":
    main()
