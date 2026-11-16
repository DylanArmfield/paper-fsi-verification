#!/usr/bin/env python3
"""SU2 inlet profile for Hron-Turek FSI1: u = 4 Umax y (H - y) / H^2, Umax = 0.3 m/s (mean 0.2).

Writes the SU2 7.5 INC velocity-inlet profile format (x, y, T, |V|, nx, ny) at every node of
the 'inlet' marker of a 2-D .su2 mesh. Same profile as the PR #75 / solids4foam codedFixedValue
(there evaluated at face centres; SU2 is vertex-based, so here at nodes).

Usage: make_inlet_profile.py MESH.su2 OUT.dat [--umax 0.3]
"""
import argparse


def inlet_nodes(mesh):
    L = open(mesh).read().split("\n")
    pts, nodes = None, set()
    i = 0
    while i < len(L):
        s = L[i]
        if s.startswith("NPOIN="):
            n = int(s.split("=")[1])
            pts = [tuple(map(float, L[i + 1 + k].split()[:2])) for k in range(n)]
            i += n
        elif s.startswith("MARKER_TAG=") and s.split("=")[1].strip() == "inlet":
            n = int(L[i + 1].split("=")[1])
            for k in range(n):
                nodes.update(map(int, L[i + 2 + k].split()[1:3]))
            i += n + 1
        i += 1
    return sorted((pts[j] for j in nodes), key=lambda p: p[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("out")
    ap.add_argument("--umax", type=float, default=0.3)
    ap.add_argument("--height", type=float, default=0.41)
    a = ap.parse_args()
    H, U = a.height, a.umax
    P = inlet_nodes(a.mesh)
    with open(a.out, "w") as f:
        f.write(f"NMARK= 1\nMARKER_TAG= inlet\nNROW={len(P)}\nNCOL=6\n")
        f.write("# COORD-X COORD-Y TEMPERATURE VELOCITY NORMAL-X NORMAL-Y\n")
        for x, y in P:
            u = max(0.0, 4.0 * U * y * (H - y) / H ** 2)
            f.write(f"{x:.15e}\t{y:.15e}\t2.881500000000000e+02\t{u:.15e}\t1.0\t0.0\n")
    # trapezoidal mean velocity on the nodes, for the record
    ys = [p[1] for p in P]
    us = [max(0.0, 4.0 * U * y * (H - y) / H ** 2) for y in ys]
    q = sum(0.5 * (us[k] + us[k + 1]) * (ys[k + 1] - ys[k]) for k in range(len(ys) - 1))
    print(f"{a.out}: {len(P)} inlet nodes, Umax {U}, trapezoidal mean {q / H:.6f} m/s (exact 2/3 Umax = {2 * U / 3:.6f})")


if __name__ == "__main__":
    main()
