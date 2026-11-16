#!/usr/bin/env python3
"""Three-level FSI1 study table (SU2 + CalculiX), from evaluate_fsi1.py --csv output.

    study_table.py RESULTS.csv            (rows for levels 1, 2, 4)

Observed order p = ln((f1 - f2)/(f2 - f4)) / ln 2 ; Richardson f_R = f4 + (f4 - f2)/(2^p - 1);
GCI(4x) = Fs |(f4 - f2)/f4| / (2^p - 1), Fs = 1.25 (as in paper-data/fsi1_four_level).
If (f1-f2)/(f2-f4) <= 0 the sequence is not monotone: no order or Richardson value is given.
The drag line is printed for completeness, but its observed order (0.20) shows it is not in the
asymptotic range, so README.md reports the 4x drag without extrapolation.
References: Featflow 7+0 (= 6+1 values used in the paper), PR #75 Richardson (OpenFOAM + CalculiX,
1x-2x-4x) and solids4foam Richardson (2x-4x-8x, paper-data/fsi1_four_level).
"""
import csv
import math
import sys

FEATFLOW = {"ux": 2.270493e-5, "uy": 8.208773e-4, "drag": 14.29426, "lift": 0.7637460}
PR75_R = {"ux": 2.2823e-5, "uy": 8.1813e-4, "drag": 14.3006, "lift": 0.76330}
S4F_R = {"ux": 2.2824e-5, "uy": 8.1825e-4, "drag": 14.3026, "lift": 0.76318}
Q = ("ux", "uy", "drag", "lift")

rows = {int(r["level"]): r for r in csv.DictReader(open(sys.argv[1]))}
assert {1, 2, 4} <= set(rows), f"need levels 1, 2, 4; have {sorted(rows)}"
f = {q: [float(rows[L][q]) for L in (1, 2, 4)] for q in Q}

print("| Quantity | 1x | 2x | 4x | Featflow 7+0 |")
print("|---|---:|---:|---:|---:|")
for q in Q:
    print(f"| {q} | " + " | ".join(f"{v:.6g}" for v in f[q]) + f" | {FEATFLOW[q]:.7g} |")
print()
print("| Quantity | Error 1x / 2x / 4x vs Featflow | Observed order | Richardson | vs Featflow | vs #75 Richardson | vs solids4foam Richardson | GCI 4x |")
print("|---|---:|---:|---:|---:|---:|---:|---:|")
for q in Q:
    f1, f2, f4 = f[q]
    errs = " / ".join(f"{100*(v/FEATFLOW[q]-1):+.2f}" for v in (f1, f2, f4)) + " %"
    r = (f1 - f2) / (f2 - f4) if f2 != f4 else float("nan")
    if r > 0 and not math.isnan(r):
        p = math.log(r) / math.log(2)
        fr = f4 + (f4 - f2) / (2 ** p - 1)
        gci = 1.25 * abs((f4 - f2) / f4) / (2 ** p - 1)
        print(f"| {q} | {errs} | {p:.2f} | {fr:.6g} | {100*(fr/FEATFLOW[q]-1):+.2f} % | "
              f"{100*(fr/PR75_R[q]-1):+.2f} % | {100*(fr/S4F_R[q]-1):+.2f} % | {100*gci:.2f} % |")
    else:
        print(f"| {q} | {errs} | not monotone (ratio {r:.2f}) | – | – | – | – | – |")
print()
print("| Level | Coupled windows | Coupling its mean / max | At cap | SU2 inner its per solve | Wall clock steady + coupled (s) |")
print("|---|---:|---:|---:|---:|---:|")
for L in (1, 2, 4):
    r = rows[L]
    print(f"| {L}x | {r['windows']} | {float(r['coupling_its_mean']):.2f} / {r['coupling_its_max']} | {r['windows_at_cap']} | "
          f"{float(r['inner_its_mean']):.1f} | {r['wall_steady_s']} + {r['wall_coupled_s']} |")
