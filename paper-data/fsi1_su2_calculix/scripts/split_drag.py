#!/usr/bin/env python3
"""Split SU2 steady drag/lift on cylinder+plate into pressure and viscous parts by integrating
surface_steady.csv (nodal Pressure and Skin_Friction_Coefficient, Cf = tau_w / (0.5 rho V_init^2) = tau_w/20)
over the marker edges of the mesh with the trapezoidal rule. Usage: split_drag.py CASEDIR [...]"""
import csv, sys
import numpy as np
from pathlib import Path

def read_mesh(p):
    L = open(p).read().split("\n"); i = 0; el = []; X = None; mk = {}
    while i < len(L):
        s = L[i]
        if s.startswith("NELEM="):
            n = int(s.split("=")[1]); el = [list(map(int, L[i+1+k].split()[1:5])) for k in range(n)]; i += n
        elif s.startswith("NPOIN="):
            n = int(s.split("=")[1]); X = np.array([list(map(float, L[i+1+k].split()[:2])) for k in range(n)]); i += n
        elif s.startswith("MARKER_TAG="):
            nm = s.split("=")[1].strip(); n = int(L[i+1].split("=")[1])
            mk[nm] = [tuple(map(int, L[i+2+k].split()[1:3])) for k in range(n)]; i += n + 1
        i += 1
    return el, X, mk

for c in map(Path, sys.argv[1:]):
    el, X, mk = read_mesh(c / "fluid-su2/fluid.su2")
    edge_el = {}
    for e in el:
        for a in range(4):
            edge_el.setdefault(frozenset((e[a], e[(a+1) % 4])), e)
    rows = list(csv.DictReader(open(c / "fluid-su2/surface_steady.csv")))
    S = {int(r["PointID"]): r for r in rows}
    tot = {}
    for m in [k for k in ("cylinder", "plate", "interface") if k in mk]:
        Fp = np.zeros(2); Fv = np.zeros(2)
        for a, b in mk[m]:
            t = X[b] - X[a]; l = np.hypot(*t)
            n = np.array([t[1], -t[0]]) / l
            cen = X[edge_el[frozenset((a, b))]].mean(0)
            if np.dot(n, cen - 0.5*(X[a]+X[b])) < 0:   # n = body outward normal, pointing INTO the fluid
                n = -n
            p = 0.5*(float(S[a]["Pressure"]) + float(S[b]["Pressure"]))
            tau = 20.0 * 0.5*np.array([float(S[a]["Skin_Friction_Coefficient_x"]) + float(S[b]["Skin_Friction_Coefficient_x"]),
                                       float(S[a]["Skin_Friction_Coefficient_y"]) + float(S[b]["Skin_Friction_Coefficient_y"])])
            Fp += -p * n * l                          # pressure traction on the body
            Fv += tau * l
        tot[m] = (Fp, Fv)
    for m, (fp, fv) in tot.items():
        print(f"  {c.name} {m:9s}: pressure drag {fp[0]:.5f} lift {fp[1]:.5f} | viscous drag {fv[0]:.5f} lift {fv[1]:.5f}")
    Fp = sum(v[0] for v in tot.values()); Fv = sum(v[1] for v in tot.values())
    print(f"{c.name}: pressure drag {Fp[0]:.5f} lift {Fp[1]:.5f} | viscous drag {Fv[0]:.5f} lift {Fv[1]:.5f} | "
          f"total drag {Fp[0]+Fv[0]:.5f} lift {Fp[1]+Fv[1]:.5f}")
