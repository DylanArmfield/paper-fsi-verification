#!/usr/bin/env python3
"""Frequency of the n = 2 ovalling vibration in ringAddedMass runs, as the solids4foam verification
driver: q(t) = (u_r(0 deg) - u_r(90 deg)) / 2 on the outer surface (CalculiX NODE PRINT of NP0 / NP90
in solid-calculix/ring.dat); zero crossings by linear interpolation, the first (t = 0) skipped; the
damped half-period is the least-squares slope of crossing time against crossing index; damping from
the peaks (parabolic interpolation), reported as the amplitude loss per period.  The wet/dry ratio
divides by a dry run on the same mesh / steps per period / alpha (--dry).

CASE is an SU2/CalculiX case folder (made by setup_case.py and run) or a published history
(histories/<run>.csv, with logs/<run>/case.json); both give the same numbers.

Usage: evaluate_ring.py CASE [CASE ...] [--dry DRYCASE] [--csv OUT]
  e.g. (from the package folder) scripts/evaluate_ring.py histories/moderate-F2.csv --dry histories/dry-F2.csv
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REF = json.loads((ROOT.parent / "febio_ring_added_mass/reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json").read_text())


def read_q(case):
    if str(case).endswith(".csv"):
        import csv
        R = list(csv.DictReader(open(case)))
        t = np.array([float(r["t"]) for r in R])
        q = np.array([(float(r["u_r_0deg"]) - float(r["u_r_90deg"])) / 2 for r in R])
        return np.concatenate([[0.0], t]), np.concatenate([[0.0], q])
    txt = (Path(case) / "solid-calculix/ring.dat").read_text()
    d = {}
    for m in re.finditer(r"displacements \(vx,vy,vz\) for set (NP0|NP90) and time\s+(\S+)\s*\n\s*\n\s*\d+\s+(\S+)\s+(\S+)", txt):
        s, t, ux, uy = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
        d.setdefault(round(t, 12), {})[s] = ux if s == "NP0" else uy
    t = np.array(sorted(k for k in d if len(d[k]) == 2))
    q = np.array([(d[k]["NP0"] - d[k]["NP90"]) / 2 for k in t])
    return np.concatenate([[0.0], t]), np.concatenate([[0.0], q])


def analyse(case):
    t, q = read_q(case)
    zc = [t[i] - q[i] * (t[i + 1] - t[i]) / (q[i + 1] - q[i])
          for i in range(1, len(q) - 1) if q[i] * q[i + 1] < 0 or (q[i + 1] == 0 and q[i] != 0)]
    k = np.arange(len(zc))
    half = np.polyfit(k, zc, 1)[0] if len(zc) > 1 else float("nan")
    omega = math.pi / half
    pk = []
    for i in range(1, len(q) - 1):
        if abs(q[i]) >= abs(q[i - 1]) and abs(q[i]) > abs(q[i + 1]):
            y0, y1, y2 = abs(q[i - 1]), abs(q[i]), abs(q[i + 1])
            den = y0 - 2 * y1 + y2
            pk.append(y1 - 0.125 * (y0 - y2) ** 2 / den if den else y1)
    loss = [1 - pk[j + 2] / pk[j] for j in range(len(pk) - 2)]
    cj = (Path(case).parent.parent / "logs" / Path(case).stem / "case.json") if str(case).endswith(".csv") else Path(case) / "case.json"
    info = json.loads(cj.read_text())
    return dict(case=str(case), level=info["level"], f=info["f"], spp=info["spp"], alpha=info["alpha"],
                omega=omega, omega_exact=info["omega_exact"], err=omega / info["omega_exact"] - 1,
                crossings=len(zc), amp=max(abs(q)), loss_per_period=float(np.mean(loss)) if loss else float("nan"),
                steps=len(t) - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--dry", default=None)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    dry = analyse(a.dry) if a.dry else None
    rows = []
    for c in a.cases:
        r = analyse(c)
        if dry and r["level"] != "dry":
            r["ratio"] = r["omega"] / dry["omega"]
            ex = REF["exact"]["levels"][r["level"]]["ratio"]
            r["ratio_err"] = r["ratio"] / ex - 1
        rows.append(r)
        s = (f"{r['case']}: {r['level']} F{r['f']} spp {r['spp']} alpha {r['alpha']}: omega {r['omega']:.7f} "
             f"(exact {r['omega_exact']:.7f}, {100 * r['err']:+.4f} %), {r['crossings']} crossings, "
             f"{r['steps']} steps, amp {r['amp']:.4e}, loss/period {100 * r['loss_per_period']:.3f} %")
        if "ratio" in r:
            s += f"; ratio {r['ratio']:.7f} ({100 * r['ratio_err']:+.4f} %)"
        print(s)
    if a.csv:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(a.csv, "w") as fh:
            fh.write(",".join(keys) + "\n")
            for r in rows:
                fh.write(",".join(str(r.get(k, "")) for k in keys) + "\n")


if __name__ == "__main__":
    main()
