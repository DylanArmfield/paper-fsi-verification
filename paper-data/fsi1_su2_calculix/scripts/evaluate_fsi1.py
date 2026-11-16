#!/usr/bin/env python3
"""Steady Hron-Turek FSI1 values of SU2 + CalculiX cases (same protocol as PR #75 evaluate.py).

    evaluate_fsi1.py CASEDIR [CASEDIR ...] [--csv OUT.csv] [--hist-dir DIR]

* Point A displacement: preCICE watch-point PointA on the solid mesh
  (solid-calculix/precice-Solid-watchpoint-PointA.log; preCICE time + 2 s = physical time).
* Drag/lift on cylinder + plate: SU2 history CD/CL, which are N/m (force reference = 1),
  converged (last) row of each time step; time window k = k-th distinct Time_Iter -> t = 2 + k*dt.
* Steady value = mean over the closing 20 % of the run (t >= 24.4 s for t_end = 30 s);
  spread = (max - min)/|mean| there; must be < 1e-3 (solids4foam / #75 test).
* Errors against Featflow 7+0 and against PR #75 (OpenFOAM + CalculiX) at the same level.
"""
import argparse
import csv
import re
from pathlib import Path

FEATFLOW = {"ux": 2.270493e-5, "uy": 8.208773e-4, "drag": 14.29426, "lift": 0.7637460}
PR75 = {  # README of PR #75, fluid and solid refined together
    1: {"ux": 2.22770e-5, "uy": 8.09851e-4, "drag": 14.21665, "lift": 0.760195},
    2: {"ux": 2.25482e-5, "uy": 8.15650e-4, "drag": 14.27168, "lift": 0.762540},
    4: {"ux": 2.26847e-5, "uy": 8.17387e-4, "drag": 14.29062, "lift": 0.763115},
}
T0 = 2.0


def read_watchpoint(p):
    lines = [l.split() for l in open(p) if l.strip()]
    head = lines[0]
    i_t, i_dx, i_dy = head.index("Time"), head.index("Displacement0"), head.index("Displacement1")
    return {round(float(r[i_t]) + T0, 9): (float(r[i_dx]), float(r[i_dy])) for r in lines[1:]}


def read_su2_history(p, dt):
    rows = list(csv.reader(open(p)))
    h = [c.strip().strip('"') for c in rows[0]]
    last = {}
    order = []
    for r in rows[1:]:
        if not r:
            continue
        ti = int(float(r[h.index("Time_Iter")]))
        if ti not in last:
            order.append(ti)
        last[ti] = r
    out, inner = {}, []
    for k, ti in enumerate(order, start=1):
        r = last[ti]
        out[round(T0 + k * dt, 9)] = (float(r[h.index("CD")]), float(r[h.index("CL")]))
    # inner iterations of every solve (all coupling iterations)
    inner = [int(float(r[h.index("Inner_Iter")])) for r in rows[1:] if r]
    plate = None
    if "CD(interface)" in h:
        r = last[order[-1]]
        plate = (float(r[h.index("CD(interface)")]), float(r[h.index("CL(interface)")]))
    return out, inner, plate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--hist-dir", default=None)
    a = ap.parse_args()
    table = []
    for c in map(Path, a.cases):
        info = (c / "CASE_INFO.txt").read_text()
        level = int(re.search(r"fluid level (\d+)", info).group(1))
        dt = float(re.search(r"dt ([0-9.e-]+) s", info).group(1))
        wp = read_watchpoint(c / "solid-calculix/precice-Solid-watchpoint-PointA.log")
        hfiles = sorted((c / "fluid-su2").glob("history_0*.csv")) or [c / "fluid-su2/history.csv"]
        fo, inner, plate = read_su2_history(hfiles[-1], dt)  # restart run: history_<RESTART_ITER>.csv
        t_end = max(wp)
        t_start = t_end - 0.2 * (t_end - T0)
        hist = {"ux": {t: v[0] for t, v in wp.items()}, "uy": {t: v[1] for t, v in wp.items()},
                "drag": {t: v[0] for t, v in fo.items()}, "lift": {t: v[1] for t, v in fo.items()}}
        its = [int(l.split()[2]) for l in open(c / "solid-calculix/precice-Solid-iterations.log")
               if l.split() and l.split()[0].isdigit()]
        wall = {k: (c / f"walltime_{k}.txt").read_text().strip() for k in ("steady", "coupled")
                if (c / f"walltime_{k}.txt").exists()}
        if "coupled" not in wall and (c / "walltime_coupled_to_last_window.txt").exists():
            wall["coupled"] = (c / "walltime_coupled_to_last_window.txt").read_text().strip() + " (to last window)"
        print(f"\n=== {c}  (level {level}, dt {dt}, last time: solid {t_end:.4f}, fluid {max(fo):.4f})")
        print(f"{'qty':5} {'mean(last 20%)':>15} {'final':>13} {'spread':>9} {'vs Featflow':>12} {'vs PR75 L'+str(level):>12}")
        row = {"case": c.name, "level": level, "dt": dt, "t_end": t_end}
        for q, h in hist.items():
            win = [v for t, v in h.items() if t >= t_start - 1e-9]
            mean = sum(win) / len(win)
            spread = (max(win) - min(win)) / abs(mean)
            e_ff = mean / FEATFLOW[q] - 1
            e_75 = mean / PR75[level][q] - 1 if level in PR75 else float("nan")
            row.update({q: mean, f"{q}_spread": spread, f"{q}_err_featflow": e_ff, f"{q}_vs_pr75": e_75})
            print(f"{q:5} {mean:15.6e} {h[max(h)]:13.6e} {spread:9.1e} {100*e_ff:11.3f}% {100*e_75:11.3f}%"
                  + ("" if spread < 1e-3 else "  NOT STEADY"))
        n_cap = sum(1 for i in its if i >= 50)
        print(f"coupled windows {len(its)}, coupling iterations mean {sum(its)/len(its):.2f}, max {max(its)}, at cap {n_cap}")
        print(f"SU2 inner iterations per solve: mean {sum(inner)/len(inner):.1f}, max {max(inner)}, solves {len(inner)}")
        if plate:
            print(f"SU2 plate-only force at end: drag {plate[0]:.6f}, lift {plate[1]:.6f} N/m")
        print(f"wall clock: {wall}")
        row.update({"windows": len(its), "coupling_its_mean": sum(its) / len(its), "coupling_its_max": max(its),
                    "windows_at_cap": n_cap, "inner_its_mean": sum(inner) / len(inner),
                    "wall_steady_s": wall.get("steady", ""), "wall_coupled_s": wall.get("coupled", "")})
        table.append(row)
        if a.hist_dir:
            Path(a.hist_dir).mkdir(parents=True, exist_ok=True)
            with open(Path(a.hist_dir) / f"{c.name}.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["t", "ux_A", "uy_A", "drag", "lift"])
                for t in sorted(set(hist["ux"]) & set(hist["drag"])):
                    w.writerow([f"{t:.6f}", hist["ux"][t], hist["uy"][t], hist["drag"][t], hist["lift"][t]])
    if a.csv and table:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(table[0]))
            w.writeheader()
            w.writerows(table)


if __name__ == "__main__":
    main()
