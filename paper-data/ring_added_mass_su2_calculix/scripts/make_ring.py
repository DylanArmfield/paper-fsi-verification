#!/usr/bin/env python3
"""Meshes and initial conditions for ringAddedMass (solids4foam tutorial definition, see
paper-data/febio_ring_added_mass/reference/s4f_ringAddedMass/verification/README.md).

Geometry: ring Ri = 0.95 .. R = 1.05 m, fluid annulus R .. b = 1.5 m.  Refinement factor F:
  solid (4 F radial) x (48 F circumferential), fluid (12 F radial) x (48 F circumferential),
  uniform in r and theta (the tutorial's cell counts); interface nodes coincide.
Initial condition (as the FEBio evidence, tools/make_ring_feb.py, V0 = 1.5e-3 m/s, a = 1):
  ring    v_r = V0 cos 2t,  v_t = V0 sin 2t (-1/2 + 3 (r - a) / (2 a))
  fluid   potential flow phi = (A r^2 + B r^-2) cos 2t with dphi/dr = V0 cos 2t at r = R and 0 at r = b
          (inviscid slip walls, so the potential field is the compatible one); p = 0 (zero ring
          acceleration at t = 0, zero displacement).

Writes:
  OUT/fluid.su2              2-D quads, markers 'interface' (r = R) and 'wall' (r = b)
  OUT/fluid_ic.csv           x, y, u, v at each fluid node (in .su2 point order)
  OUT/mesh.inp               C3D8I, one layer of depth 1 m (u_z = 0 -> plane strain); sets Nall,
                             Nsurface (outer surface, the interface), NP0 / NP90 (outer nodes at 0 / 90 deg)
  OUT/initial_velocity.inp   *INITIAL CONDITIONS, TYPE=VELOCITY for every solid node

Usage: make_ring.py F OUTDIR [--v0 1.5e-3]
"""
import argparse
import math
from pathlib import Path

RI, RO, B, A = 0.95, 1.05, 1.5, 1.0


def potential_coeffs(v0):
    # dphi/dr = 2 A r - 2 Bc r^-3 :  = v0 at RO, 0 at B  ->  Bc = A B^4
    a = v0 / (2.0 * (RO - B ** 4 * RO ** -3))
    return a, a * B ** 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("f", type=int)
    ap.add_argument("out")
    ap.add_argument("--v0", type=float, default=1.5e-3)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    nt, nrs, nrf = 48 * a.f, 4 * a.f, 12 * a.f
    th = [2 * math.pi * i / nt for i in range(nt)]

    # ---------------- fluid ----------------
    rf = [RO + (B - RO) * j / nrf for j in range(nrf + 1)]
    fid = lambda i, j: (i % nt) + nt * j
    pts, ic = [], []
    Ac, Bc = potential_coeffs(a.v0)
    for j in range(nrf + 1):
        for i in range(nt):
            r, t = rf[j], th[i]
            x, y = r * math.cos(t), r * math.sin(t)
            vr = 2 * (Ac * r - Bc / r ** 3) * math.cos(2 * t)
            vt = -2 * (Ac * r + Bc / r ** 3) * math.sin(2 * t)
            u = vr * math.cos(t) - vt * math.sin(t)
            v = vr * math.sin(t) + vt * math.cos(t)
            pts.append((x, y))
            ic.append((x, y, u, v))
    # quads counter-clockwise: (i, j) -> (i, j+1) -> (i+1, j+1) -> (i+1, j)  [r outward, theta CCW]
    els = [(fid(i, j), fid(i, j + 1), fid(i + 1, j + 1), fid(i + 1, j)) for j in range(nrf) for i in range(nt)]
    for e in els[:1]:
        (x0, y0), (x1, y1), (x2, y2) = pts[e[0]], pts[e[1]], pts[e[2]]
        assert (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0) > 0, "quad orientation"
    L = ["NDIME= 2", f"NELEM= {len(els)}"]
    L += [f"9 {e[0]} {e[1]} {e[2]} {e[3]} {k}" for k, e in enumerate(els)]
    L.append(f"NPOIN= {len(pts)}")
    L += [f"{x:.16e} {y:.16e} {k}" for k, (x, y) in enumerate(pts)]
    L.append("NMARK= 2")
    for name, j in (("interface", 0), ("wall", nrf)):
        L += [f"MARKER_TAG= {name}", f"MARKER_ELEMS= {nt}"]
        L += [f"3 {fid(i, j)} {fid(i + 1, j)}" for i in range(nt)]
    (out / "fluid.su2").write_text("\n".join(L) + "\n")
    (out / "fluid_ic.csv").write_text("x,y,u,v\n" + "".join(f"{x:.16e},{y:.16e},{u:.16e},{v:.16e}\n" for x, y, u, v in ic))

    # ---------------- solid ----------------
    rs = [RI + (RO - RI) * j / nrs for j in range(nrs + 1)]
    sid = lambda i, j, k: 1 + (i % nt) + nt * (j + (nrs + 1) * k)
    L = ["*NODE, NSET=Nall"]
    V = ["*INITIAL CONDITIONS, TYPE=VELOCITY"]
    surf = []
    for k in range(2):
        for j in range(nrs + 1):
            for i in range(nt):
                r, t = rs[j], th[i]
                n = sid(i, j, k)
                L.append(f"{n}, {r * math.cos(t):.12e}, {r * math.sin(t):.12e}, {float(k):.1f}")
                vr = a.v0 * math.cos(2 * t)
                vt = a.v0 * math.sin(2 * t) * (-0.5 + 1.5 * (r - A) / A)
                V.append(f"{n}, 1, {vr * math.cos(t) - vt * math.sin(t):.12e}")
                V.append(f"{n}, 2, {vr * math.sin(t) + vt * math.cos(t):.12e}")
                if j == nrs:
                    surf.append(n)
    L.append("*ELEMENT, TYPE=C3D8I, ELSET=Eall")
    e = 0
    for j in range(nrs):
        for i in range(nt):
            e += 1
            q = [(i, j), (i, j + 1), (i + 1, j + 1), (i + 1, j)]
            L.append(f"{e}, " + ", ".join(str(sid(p, s, k)) for k in (0, 1) for p, s in q))
    for name, ids in (("Nsurface", surf), ("NP0", [sid(0, nrs, 0), sid(0, nrs, 1)]),
                      ("NP90", [sid(nt // 4, nrs, 0), sid(nt // 4, nrs, 1)])):
        L.append(f"*NSET, NSET={name}")
        L += [", ".join(map(str, ids[b:b + 10])) + "," for b in range(0, len(ids), 10)]
    (out / "mesh.inp").write_text("\n".join(L) + "\n")
    (out / "initial_velocity.inp").write_text("\n".join(V) + "\n")
    print(f"F {a.f}: fluid {nrf}x{nt} = {len(els)} quads, {len(pts)} nodes; solid {nrs}x{nt} = {e} C3D8I; "
          f"{nt} interface points; V0 {a.v0}; potential A {Ac:.6e} B {Bc:.6e}")


if __name__ == "__main__":
    main()
