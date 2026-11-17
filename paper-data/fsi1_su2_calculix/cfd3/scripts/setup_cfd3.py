#!/usr/bin/env python3
"""Create one SU2 Hron-Turek CFD3 case (rigid flag, unsteady) for the beta / cost study.

    setup_cfd3.py LEVEL MESH.su2 OUTDIR [--dt DT] [--t-end 10] [--beta 4.1] [--inner-iter 1000]
                  [--inner-minval -9] [--write-every 500] [--restart-from CASEDIR --restart-time T]
                  [--set KEY=VALUE ...]

MESH.su2: the standard TurekHron2006 family x LEVEL made by ../../scripts/make_fluid_mesh.py from
the #75 blockMeshDict (the FSI1 meshes: 5336 / 21344 / 85376 cells), with the plate marker kept as
'plate' (no coupling). dt default: the paper's matched CFD3 sequence (1e-3, 5e-4, 2.5e-4).
Inlet: parabolic, Umax 3 m/s (mean 2), at the inlet nodes, impulsive.

--restart-from copies the two BDF2 restart files of CASEDIR at time T (restart_flow_NNNNN.dat for
the time iterations N-1 and N) and continues from there; TIME_ITER is then the absolute end iteration.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # cfd3/
FSI1_SCRIPTS = ROOT.parent / "scripts"  # make_inlet_profile.py, shared with FSI1
DT = {1: 1e-3, 2: 5e-4, 4: 2.5e-4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("level", type=int, choices=sorted(DT))
    ap.add_argument("mesh")
    ap.add_argument("outdir")
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--t-end", type=float, default=10.0)
    ap.add_argument("--beta", type=float, default=4.1)
    ap.add_argument("--inner-iter", type=int, default=1000)
    ap.add_argument("--inner-minval", type=float, default=-9.0)
    ap.add_argument("--write-every", type=int, default=500)
    ap.add_argument("--restart-from", default=None)
    ap.add_argument("--restart-time", type=float, default=None)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override/add a cfg option (solver settings for diagnostics), e.g. LINEAR_SOLVER_ITER=200")
    a = ap.parse_args()
    L = a.level
    dt = a.dt or DT[L]
    out = Path(a.outdir).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    out.mkdir(parents=True)
    mesh = Path(a.mesh).resolve()
    m = mesh.read_text()
    assert m.count("MARKER_TAG= plate\n") == 1
    (out / "fluid.su2").write_text(m)
    subprocess.run([sys.executable, "-I", FSI1_SCRIPTS / "make_inlet_profile.py", out / "fluid.su2",
                    out / "inlet_profile.dat", "--umax", "3.0"], check=True)

    n_end = round(a.t_end / dt)
    restart, rit, note = "NO", 0, ""
    if a.restart_from:
        src = Path(a.restart_from).resolve()
        n0 = round(a.restart_time / dt)
        for n in (n0 - 1, n0):
            f = src / f"restart_flow_{n:05d}.dat"
            if not f.exists():
                sys.exit(f"missing {f}")
            shutil.copy2(f, out / f.name)
        restart, rit = "YES", n0 + 1
        note = f"; restart from {src.name} at t = {a.restart_time} (iter {n0})"
    # unsteady SU2 reads the inlet profile with a time-iteration suffix (as ../../scripts/run-case.sh);
    # a restart at iteration R reads inlet_profile_<R-2>.dat
    for n in sorted({0, rit, max(0, rit - 1), max(0, rit - 2)}):
        shutil.copy2(out / "inlet_profile.dat", out / f"inlet_profile_{n:05d}.dat")
    cfg = (ROOT / "case/cfd3.cfg").read_text()
    sub = dict(RESTART=restart, RESTART_ITER=rit, DT=repr(dt), TIME_ITER=n_end, INNER_ITER=a.inner_iter,
               BETA=repr(a.beta), INNER_MINVAL=a.inner_minval, WRT=a.write_every)
    for k, v in sub.items():
        cfg = cfg.replace(f"@{k}@", str(v))
    for kv in a.set:
        k, v = kv.split("=", 1)
        k, v = k.strip(), v.strip()
        cfg, n = re.subn(rf"(?m)^{k}=.*$", f"{k}= {v}", cfg)
        if n == 0:
            cfg += f"{k}= {v}\n"
        note += f"; {k}= {v}"
    if re.search(r"@[A-Z_]+@", cfg):
        sys.exit("unfilled placeholder in cfd3.cfg")
    (out / "cfd3.cfg").write_text(cfg)
    ncell = int(re.search(r"NELEM=\s*(\d+)", m).group(1))
    info = (f"CFD3 level {L} ({mesh.name}, {ncell} cells), dt {dt} s, t_end {a.t_end} s ({n_end} steps), "
            f"BETA_FACTOR {a.beta}, inner iter cap {a.inner_iter}, inner minval {a.inner_minval}{note}")
    (out / "CASE_INFO.txt").write_text(info + "\n")
    print(info + f"\ncase in {out}")


if __name__ == "__main__":
    main()
