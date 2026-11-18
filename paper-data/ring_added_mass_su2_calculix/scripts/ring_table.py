#!/usr/bin/env python3
"""Summary tables of the SU2 + CalculiX ringAddedMass runs, from the published histories and logs only.

Writes (in the package folder):
  ring_su2_calculix.csv   one row per run: settings (logs/<run>/CASE_INFO.txt, case.json), cores, wall time,
                          exit status, frequency and its error, wet/dry ratio and its error (dry run on the same
                          mesh and steps per period), amplitude loss per period, coupling iterations per window,
                          SU2 inner iterations and the median final log10 rms residuals of the SU2 solves
  ring_mesh_study.csv     ratio error at F = 1, 2, 4 (100 steps per period), observed order, Richardson ratio
  ring_time_study.csv     dry and moderate level, F = 2, 25 / 50 / 100 steps per period
  ring_acceptance.csv     the solids4foam verification acceptance criteria (reference JSON tolerances) as
                          they apply to this solid (see README)
  beta/beta_single_step.csv  BETA_FACTOR single-step check (beta/logs/*.log)

Usage (from the package folder):  python3 scripts/ring_table.py
"""
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import evaluate_ring as ev  # noqa: E402

TOL = ev.REF["tolerances"]
RUNS = ["dry-F1", "dry-F2", "dry-F4", "dry-F2-spp25", "dry-F2-spp50",
        "weak-F1", "weak-F2", "weak-F4", "moderate-F1", "moderate-F2", "moderate-F4",
        "moderate-F2-spp25", "moderate-F2-spp50", "moderate-F1-beta1e5", "moderate-F1-beta1e8",
        "strong-F1", "strong-F2", "strong-F4", "strong-F1-reuse10"]
LEVELS = ["weak", "moderate", "strong"]
COUPLING_CAP = 50   # preCICE max-iterations
BETA_LOGS = {"n_50.log": (4.1, 50), "n_300.log": (4.1, 300), "n_2000.log": (4.1, 2000),
             "v_b1e4.log": (1e4, 300), "v_b1e6.log": (1e6, 300),
             "v_b1E5.log": (1e5, 1000), "v_b1E6.log": (1e6, 1000), "v_b1E7.log": (1e7, 1000),
             "v_b1E8.log": (1e8, 1000), "v_b1E10.log": (1e10, 1000)}


def fmt(v):
    return f"{v:.7g}" if isinstance(v, float) else ("" if v is None else v)


def write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([fmt(r.get(k)) for k in keys])
    print(f"wrote {path.relative_to(ROOT)}")


def dry_of(cj):
    return f"dry-F{cj['f']}" + ("" if cj["spp"] == 100 else f"-spp{cj['spp']}")


def run_row(run):
    L = ROOT / "logs" / run
    cj = json.loads((L / "case.json").read_text())
    info = (L / "CASE_INFO.txt").read_text()
    head = (L / "slurm.out").read_text().splitlines()[0]
    r = dict(run=run, level=cj["level"], F=cj["f"], steps_per_period=cj["spp"], dt=cj["dt"], hht_alpha=cj["alpha"])
    if cj["level"] != "dry":
        r["beta_factor"] = float(re.search(r"BETA_FACTOR ([0-9.eE+-]+)", info).group(1))
        r["iqn_reuse"] = int(re.search(r'time-windows-reused value="(\d+)"', (L / "precice-config.xml").read_text()).group(1))
        r["coupling_rel_tol"] = cj["ctol"]
        r["su2_inner_cap"] = int(re.search(r"inner cap (\d+)", info).group(1))
    r.update(slurm_job=int(re.search(r"Job (\d+)", head).group(1)), cores=int(re.search(r"tasks (\d+)", head).group(1)),
             exit_status=int((L / "exit_status.txt").read_text()), wall_s=int((L / "walltime.txt").read_text()))
    r["core_hours"] = r["wall_s"] * r["cores"] / 3600
    hist = ROOT / "histories" / f"{run}.csv"
    H = list(csv.DictReader(open(hist)))
    r["steps_run"] = len(H)
    r["steps_planned"] = cj["nsteps"]
    if len(H) == cj["nsteps"]:
        a = ev.analyse(str(hist))
        r.update(omega=a["omega"], omega_exact=a["omega_exact"], omega_err_pct=100 * a["err"], zero_crossings=a["crossings"],
                 amplitude=a["amp"], loss_per_period_pct=100 * a["loss_per_period"])
        if cj["level"] != "dry":
            d = ev.analyse(str(ROOT / "histories" / f"{dry_of(cj)}.csv"))
            ex = ev.REF["exact"]["levels"][cj["level"]]["ratio"]
            r.update(dry_run=dry_of(cj), ratio=a["omega"] / d["omega"], ratio_exact=ex,
                     ratio_err_pct=100 * (a["omega"] / d["omega"] / ex - 1))
    if cj["level"] != "dry":
        its = [int(h["coupling_iters"]) for h in H if h["coupling_iters"]]
        r.update(coupling_iters_mean=statistics.mean(its), coupling_iters_max=max(its),
                 windows_at_coupling_cap=sum(i >= COUPLING_CAP for i in its))
        S = list(csv.DictReader(open(ROOT / "histories" / f"{run}_su2.csv")))
        r.update(su2_solves=len(S), su2_solves_at_inner_cap=sum(int(s["inner_iters"]) >= r["su2_inner_cap"] for s in S))
        for c in ("p", "u", "v"):
            r[f"su2_median_final_log10_rms_{c}"] = statistics.median(float(s[f"log10_rms_{c}"]) for s in S)
    return r


def main():
    rows = {run: run_row(run) for run in RUNS}
    write(ROOT / "ring_su2_calculix.csv", list(rows.values()))

    mesh = []
    for lev in LEVELS:
        e = [rows[f"{lev}-F{f}"]["ratio_err_pct"] for f in (1, 2, 4)]
        rt = [rows[f"{lev}-F{f}"]["ratio"] for f in (1, 2, 4)]
        p = math.log2((e[0] - e[1]) / (e[1] - e[2]))
        rr = rt[2] + (rt[2] - rt[1]) / (2 ** p - 1)
        ex = ev.REF["exact"]["levels"][lev]["ratio"]
        mesh.append(dict(level=lev, F1_ratio_err_pct=e[0], F2_ratio_err_pct=e[1], F4_ratio_err_pct=e[2], observed_order=p,
                         ratio_richardson=rr, ratio_exact=ex, richardson_ratio_err_pct=100 * (rr / ex - 1),
                         F4_omega_err_pct=rows[f"{lev}-F4"]["omega_err_pct"],
                         F2_loss_per_period_pct=rows[f"{lev}-F2"]["loss_per_period_pct"],
                         F4_loss_per_period_pct=rows[f"{lev}-F4"]["loss_per_period_pct"],
                         coupling_iters_mean_F4=rows[f"{lev}-F4"]["coupling_iters_mean"]))
    write(ROOT / "ring_mesh_study.csv", mesh)

    time = []
    for lev in ("dry", "moderate"):
        R = [rows[f"{lev}-F2-spp{s}"] for s in (25, 50)] + [rows[f"{lev}-F2"]]
        w = [x["omega"] for x in R]
        q = math.log2((w[1] - w[0]) / (w[2] - w[1]))
        wr = w[2] + (w[2] - w[1]) / 3   # order-2 Richardson of 50 and 100 steps per period (as solids4foam)
        t = dict(level=lev, F=2, omega_err_pct_spp25=R[0]["omega_err_pct"], omega_err_pct_spp50=R[1]["omega_err_pct"],
                 omega_err_pct_spp100=R[2]["omega_err_pct"], observed_time_order=q,
                 omega_time_richardson_err_pct=100 * (wr / R[2]["omega_exact"] - 1))
        if lev != "dry":
            t.update(ratio_err_pct_spp25=R[0]["ratio_err_pct"], ratio_err_pct_spp50=R[1]["ratio_err_pct"],
                     ratio_err_pct_spp100=R[2]["ratio_err_pct"],
                     loss_per_period_pct_spp25=R[0]["loss_per_period_pct"], loss_per_period_pct_spp50=R[1]["loss_per_period_pct"],
                     loss_per_period_pct_spp100=R[2]["loss_per_period_pct"])
        time.append(t)
    write(ROOT / "ring_time_study.csv", time)

    acc = []
    def chk(name, value, limit, ok, applies="yes", note=""):
        acc.append(dict(check=name, value=value, limit=limit, passed=("yes" if ok else "no") if applies == "yes" else "",
                        applies=applies, note=note))
    for m in mesh:
        chk(f"ratio error, finest mesh (F4), {m['level']}", m["F4_ratio_err_pct"], 100 * TOL["ratio"],
            abs(m["F4_ratio_err_pct"]) <= 100 * TOL["ratio"])
        dec = m["F1_ratio_err_pct"] > m["F2_ratio_err_pct"] > m["F4_ratio_err_pct"] > 0
        chk(f"ratio error decreases under refinement, {m['level']}", "", "", dec)
        chk(f"observed mesh order of the ratio, {m['level']}", m["observed_order"], TOL["minMeshOrder"],
            m["observed_order"] >= TOL["minMeshOrder"], note="solids4foam: standard-solid criterion")
    for lev in ["dry"] + LEVELS:
        v = rows[f"{lev}-F4"]["omega_err_pct"]
        chk(f"frequency error, finest mesh (F4), 100 steps/period, {lev}", v, 100 * TOL["omega"], abs(v) <= 100 * TOL["omega"])
    for lev in LEVELS:
        v = rows[f"{lev}-F2"]["loss_per_period_pct"]
        chk(f"amplitude loss per period, F2, 100 steps/period, {lev}", v, 100 * TOL["maxLossPerPeriod"],
            v <= 100 * TOL["maxLossPerPeriod"])
    gmax = max(-r["loss_per_period_pct"] for r in rows.values() if "loss_per_period_pct" in r)
    chk("largest amplitude growth per period, any run", gmax, 100 * TOL["maxGrowthPerPeriod"], gmax <= 100 * TOL["maxGrowthPerPeriod"])
    for t in time:
        chk(f"observed time order of the frequency, F2, {t['level']}", t["observed_time_order"], TOL["minTimeOrder"],
            t["observed_time_order"] >= TOL["minTimeOrder"])
        chk(f"time-extrapolated frequency error, F2, {t['level']}", t["omega_time_richardson_err_pct"], 100 * TOL["omegaExtrapolated"],
            False, applies="no", note="contains the F2 spatial error of the C3D8I solid; the solids4foam tolerance assumes "
            "a solid converged in space on every mesh (its high-order solid)")
    chk("frequency change between meshes", "", 100 * TOL["meshChange"], False, applies="no",
        note="high-order-solid criterion; the C3D8I solid is not converged in space on F1 (dry F1 error in ring_su2_calculix.csv)")
    chk("Robin-Neumann and IQN-ILS frequencies agree", "", 100 * TOL["couplingAgreement"], False, applies="no",
        note="only IQN-ILS was run")
    write(ROOT / "ring_acceptance.csv", acc)

    beta = []
    for f, (b, cap) in BETA_LOGS.items():
        txt = (ROOT / "beta/logs" / f).read_text()
        last = [ln for ln in txt.splitlines() if re.match(r"^\|\s+2\|", ln)][-1].strip("|").split("|")
        F2, F0 = map(float, re.search(r"MODAL F2 (\S+)\s+F0\(uniform\) (\S+)", txt).groups())
        beta.append(dict(log=f, beta_factor=b, inner_cap=cap, inner_iters=int(last[2]) + 1, final_log10_rms_p=float(last[3]),
                         final_log10_rms_u=float(last[4]), final_log10_rms_v=float(last[5]), modal_force_F2=F2, uniform_force_F0=F0))
    ref = next(x["modal_force_F2"] for x in beta if x["log"] == "v_b1E6.log")
    for x in beta:
        x["F2_over_beta1e6"] = x["modal_force_F2"] / ref
    write(ROOT / "beta/beta_single_step.csv", beta)

    for r in rows.values():
        s = f"{r['run']:20s} job {r['slurm_job']} {r['cores']} cores {r['wall_s']:6d} s exit {r['exit_status']}"
        if "omega" in r:
            s += f": omega {r['omega']:.6f} ({r['omega_err_pct']:+.3f} %), loss {r['loss_per_period_pct']:.3f} %/period"
        if "ratio" in r:
            s += f", ratio {r['ratio_err_pct']:+.4f} %"
        if "coupling_iters_mean" in r:
            s += (f"; coupling {r['coupling_iters_mean']:.2f} (max {r['coupling_iters_max']}); SU2 solves at cap "
                  f"{r['su2_solves_at_inner_cap']}/{r['su2_solves']}, median rms p/u {r['su2_median_final_log10_rms_p']:.2f}"
                  f"/{r['su2_median_final_log10_rms_u']:.2f}")
        if r["steps_run"] < r["steps_planned"]:
            s += f"; stopped after {r['steps_run']} of {r['steps_planned']} steps"
        print(s)
    for m in mesh:
        print(f"{m['level']:9s} ratio err F1/F2/F4 {m['F1_ratio_err_pct']:+.4f} / {m['F2_ratio_err_pct']:+.4f} / "
              f"{m['F4_ratio_err_pct']:+.4f} %, order {m['observed_order']:.2f}, Richardson {m['richardson_ratio_err_pct']:+.4f} %")
    for a in acc:
        print(f"[{a['passed'] or 'n/a':3s}] {a['check']}: {fmt(a['value'])} (limit {a['limit']}) {a['note']}")
    for x in beta:
        print(f"beta {x['beta_factor']:g} cap {x['inner_cap']}: rms p {x['final_log10_rms_p']:.2f}, F2 {x['modal_force_F2']:.4f} "
              f"({x['F2_over_beta1e6']:.3f} of beta 1e6)")


if __name__ == "__main__":
    main()
