#!/usr/bin/env python3
"""Create one Hron-Turek FSI1 preCICE case (SU2 incompressible + CalculiX) for a mesh level.

    setup_case.py LEVEL OUTDIR [--solid-level S] [--depth 1.0]
                  [--inner-iter 200] [--inner-minval -9] [--steady-iter 20000]

Mirrors PR #75 (paper-data/fsi1_precice_calculix) scripts/setup-case.py so that only the fluid code differs:
  * fluid mesh: the TurekHron2006 blockMesh x LEVEL (make_fluid_mesh.py, node-identical to #75);
  * solid: (105*S) x (6*S) x 1 C3D8I, x 0.2489898-0.6, y 0.19-0.21, Nfix at x0, Nsurface =
    top+bottom+tip nodes -- generator copied from #75; slab depth = --depth (default 1 m,
    because SU2 forces are per metre depth and the quasi-2D CalculiX adapter applies the 2-D
    nodal force to the whole slab; #75 uses 0.015 m with OpenFOAM's 0.015 m-deep forces);
  * dt = 0.025 s (L1, L2), 0.0125 s (L4); coupled 28 s.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent / "case"
# the PR #75 fluid blockMesh (TurekHron2006), already in this repository
BLOCKMESH = HERE.parent.parent / "fsi1_precice_calculix/case/fluid-openfoam/system/blockMeshDict"
DT = {1: 0.025, 2: 0.025, 4: 0.0125}


def solid_mesh(slevel, depth):
    # PR #75 generator, depth as parameter
    x0, x1, y0, y1 = 0.2489898, 0.6, 0.19, 0.21
    nx, ny = 105 * slevel, 6 * slevel
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    lines = ["*NODE, NSET=Nall"]
    fix, surf = [], []
    for k in range(2):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x = x0 + (x1 - x0) * i / nx
                y = y0 + (y1 - y0) * j / ny
                lines.append(f"{nid(i,j,k)}, {x:.12e}, {y:.12e}, {depth*k:.12e}")
                if i == 0:
                    fix.append(nid(i, j, k))
                if j == 0 or j == ny or i == nx:
                    surf.append(nid(i, j, k))
    lines.append("*ELEMENT, TYPE=C3D8I, ELSET=Eall")
    e = 0
    for j in range(ny):
        for i in range(nx):
            e += 1
            n = [nid(i, j, 0), nid(i+1, j, 0), nid(i+1, j+1, 0), nid(i, j+1, 0),
                 nid(i, j, 1), nid(i+1, j, 1), nid(i+1, j+1, 1), nid(i, j+1, 1)]
            lines.append(f"{e}, " + ", ".join(map(str, n)))

    def nset(name, ids):
        lines.append(f"*NSET, NSET={name}")
        for a in range(0, len(ids), 10):
            lines.append(", ".join(map(str, ids[a:a + 10])) + ",")
    nset("Nfix", fix)
    nset("Nsurface", surf)
    return "\n".join(lines) + "\n", nx, ny, e, len(surf) // 2


def subst(path, **kw):
    t = path.read_text()
    for k, v in kw.items():
        t = t.replace(f"@{k}@", str(v))
    path.write_text(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("level", type=int, choices=sorted(DT))
    ap.add_argument("outdir")
    ap.add_argument("--solid-level", type=int, default=None)
    ap.add_argument("--depth", type=float, default=1.0)
    ap.add_argument("--inner-iter", type=int, default=200)
    ap.add_argument("--inner-minval", type=float, default=-9.0)  # -7 left a ~1e-6 force-noise floor vs the 1e-6 coupling tol
    ap.add_argument("--steady-iter", type=int, default=20000)
    a = ap.parse_args()
    L, S = a.level, a.solid_level or a.level
    out = Path(a.outdir).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    shutil.copytree(TEMPLATE, out)
    dt = DT[L]

    # fluid
    fl = out / "fluid-su2"
    mesh = out / f"fsi1_L{L}.su2"   # generated with OpenFOAM blockMesh in out/foam-L<L>
    subprocess.run([sys.executable, "-I", HERE / "make_fluid_mesh.py", BLOCKMESH, str(L), mesh], check=True)
    # the SU2 preCICE adapter hard-codes the coupled marker name "interface": rename "plate"
    m = mesh.read_text()
    assert m.count("MARKER_TAG= plate\n") == 1
    (fl / "fluid.su2").write_text(m.replace("MARKER_TAG= plate\n", "MARKER_TAG= interface\n"))
    subprocess.run([sys.executable, "-I", HERE / "make_inlet_profile.py", fl / "fluid.su2",
                    fl / "inlet_profile.dat"], check=True)
    subst(fl / "steady_rigid.cfg", MESH="fluid.su2", STEADY_ITER=a.steady_iter)
    subst(fl / "coupled.cfg", MESH="fluid.su2", DT=repr(dt), INNER_ITER=a.inner_iter,
          INNER_MINVAL=a.inner_minval)

    # solid
    sd = out / "solid-calculix"
    txt, nx, ny, ne, nif = solid_mesh(S, a.depth)
    (sd / "mesh.inp").write_text(txt)
    subst(sd / "plate.inp", DT=repr(dt))

    # coupling
    subst(out / "precice-config.base.xml", DT=repr(dt))

    for f in out.rglob("*"):
        if f.is_file() and f.suffix in (".cfg", ".inp", ".xml") and "@" in f.read_text():
            sys.exit(f"unfilled placeholder in {f}")
    info = (f"fluid level {L} ({mesh.name}), solid level {S}: {nx}x{ny} C3D8I ({ne} elements, "
            f"{nif} interface points), depth {a.depth} m, dt {dt} s, inner iter {a.inner_iter}, "
            f"inner minval {a.inner_minval}, steady iter {a.steady_iter}")
    (out / "CASE_INFO.txt").write_text(info + "\n")
    print(info + f"\ncase in {out}")


if __name__ == "__main__":
    main()
