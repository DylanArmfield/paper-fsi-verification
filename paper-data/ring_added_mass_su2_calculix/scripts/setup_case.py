#!/usr/bin/env python3
"""Create one ringAddedMass case: dry (CalculiX alone) or wet (SU2 + CalculiX via preCICE).

    setup_case.py LEVEL F OUTDIR [--spp 100] [--periods 4] [--alpha -0.05] [--ctol 1e-6] [--reuse 10]
                  [--inner-iter 500] [--inner-minval -12] [--beta 4.1] [--write-every 25]

LEVEL dry | weak | moderate | strong (rho_f 14 / 140 / 1400); F mesh refinement factor (1, 2, 4).
dt = T / spp with T = 2 pi / omega_exact of the level (as the solids4foam verification driver:
fixed steps per period, so the time error is similar at every level), run for --periods periods.
Exact frequencies from the solids4foam verification reference JSON (copied in ../febio_ring_added_mass/reference/).
The published runs: see ../README.md (Reproduction).
"""
import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF = json.loads((ROOT.parent / "febio_ring_added_mass/reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json").read_text())
RHOF = {"weak": 14.0, "moderate": 140.0, "strong": 1400.0}


def subst(path, **kw):
    t = path.read_text()
    for k, v in kw.items():
        t = t.replace(f"@{k}@", str(v))
    path.write_text(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("level", choices=["dry", "weak", "moderate", "strong"])
    ap.add_argument("f", type=int)
    ap.add_argument("outdir")
    ap.add_argument("--spp", type=int, default=100)
    ap.add_argument("--periods", type=float, default=4)
    ap.add_argument("--alpha", type=float, default=-0.05)
    ap.add_argument("--ctol", type=float, default=1e-6)
    ap.add_argument("--reuse", type=int, default=10, help="IQN-ILS time-windows-reused")
    ap.add_argument("--inner-iter", type=int, default=500)
    ap.add_argument("--inner-minval", type=float, default=-12)
    ap.add_argument("--beta", type=float, default=4.1)
    ap.add_argument("--write-every", type=int, default=25)
    ap.add_argument("--screen-inner", type=int, default=1000)
    ap.add_argument("--v0", type=float, default=1.5e-3)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override/add an SU2 cfg option")
    a = ap.parse_args()
    out = Path(a.outdir).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    omega = REF["exact"]["dry"]["omega"] if a.level == "dry" else REF["exact"]["levels"][a.level]["omega"]
    T = 2 * math.pi / omega
    dt = T / a.spp
    nsteps = round(a.periods * a.spp)
    tend = nsteps * dt

    shutil.copytree(ROOT / "case", out)
    sd, fl = out / "solid-calculix", out / "fluid-su2"
    r = subprocess.run([sys.executable, "-I", HERE / "make_ring.py", str(a.f), str(sd), "--v0", repr(a.v0)],
                       check=True, capture_output=True, text=True)
    subst(sd / "ring.inp", ALPHA=repr(a.alpha), DT=repr(dt), TEND=repr(tend), WRT=a.write_every)
    if a.level == "dry":
        shutil.rmtree(fl)
        (out / "precice-config.base.xml").unlink()
        (sd / "config.yml").unlink()
        for f in ("fluid.su2", "fluid_ic.csv"):
            (sd / f).unlink()
    else:
        for f in ("fluid.su2", "fluid_ic.csv"):
            shutil.move(sd / f, fl / f)
        # potential-flow initial condition as an ASCII SU2 restart (format of SU2 7.5.1 INC, no energy eq.)
        rows = list(csv.DictReader(open(fl / "fluid_ic.csv")))
        hdr = '"PointID","x","y","Pressure","Velocity_x","Velocity_y"\n'
        body = "".join(f"{k}, {float(x['x']):.15e}, {float(x['y']):.15e}, 0.000000000000000e+00, "
                       f"{float(x['u']):.15e}, {float(x['v']):.15e}\n" for k, x in enumerate(rows))
        (fl / "restart_ic.csv").write_text(hdr + body)
        subst(fl / "coupled.cfg", RHOF=RHOF[a.level], BETA=repr(a.beta), DT=repr(dt), INNER_ITER=a.inner_iter,
              INNER_MINVAL=a.inner_minval, WRT=a.write_every, SCREEN_INNER=a.screen_inner)
        for kv in a.set:
            k, v = (x.strip() for x in kv.split("=", 1))
            cfg = (fl / "coupled.cfg").read_text()
            cfg, n = re.subn(rf"(?m)^{k}=.*$", f"{k}= {v}", cfg)
            (fl / "coupled.cfg").write_text(cfg if n else cfg + f"{k}= {v}\n")
        subst(out / "precice-config.base.xml", DT=repr(dt), NSTEPS=nsteps, WRT=a.write_every, CTOL=repr(a.ctol), REUSE=a.reuse)
    for f in out.rglob("*"):
        if f.is_file() and f.suffix in (".cfg", ".inp", ".xml", ".yml") and re.search(r"@[A-Z_]+@", f.read_text()):
            sys.exit(f"unfilled placeholder in {f}")
    info = (f"ringAddedMass {a.level} (rho_f {RHOF.get(a.level, '-')}), mesh factor {a.f}: {r.stdout.strip()}; "
            f"omega_exact {omega:.9f} rad/s, T {T:.9f} s, {a.spp} steps/period, dt {dt:.9e} s, {nsteps} steps, "
            f"t_end {tend:.9f} s; HHT alpha {a.alpha}")
    if a.level != "dry":
        info += (f"; coupling rel tol {a.ctol}; SU2 inner cap {a.inner_iter}, minval {a.inner_minval}, "
                 f"BETA_FACTOR {a.beta}" + (f"; IQN-ILS time-windows-reused {a.reuse}" if a.reuse != 10 else "") + "".join(f"; {x}" for x in a.set))
    (out / "CASE_INFO.txt").write_text(info + "\n")
    (out / "case.json").write_text(json.dumps(dict(level=a.level, f=a.f, omega_exact=omega, T=T, dt=dt,
                                                   nsteps=nsteps, spp=a.spp, alpha=a.alpha, ctol=a.ctol), indent=1))
    print(info + f"\ncase in {out}")


if __name__ == "__main__":
    main()
