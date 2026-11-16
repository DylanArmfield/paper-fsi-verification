#!/usr/bin/env python3
"""SU2 (pysu2) with a PRESCRIBED plate motion, no preCICE: isolates the SU2 adapter's
checkpoint (SaveOldState / ReloadOldState) from the coupling.

Plate displacement: u_y(x) = a_n * ((x - x0)/(x1 - x0))^2, a_n = AMP * n (n = time step 1, 2, ...).
Modes (same sequence of calls as run/SU2_preCICE_FSI.py):
  single      : per step  set disp(a_n) -> Preprocess/Run/Postprocess/Update/Monitor -> Output
  reload-same : per step  Save; solve with disp(a_n); Reload; solve with disp(a_n) again; Output
                (an exact checkpoint must reproduce 'single' to round-off)
  reload-iqn  : per step  Save; solve with disp(0.5 a_n); Reload; solve with disp(a_n); Output
                (mimics an under-relaxed first coupling iteration)

Usage: mpirun -np N python3 prescribed.py CFG MODE NSTEPS [AMP]
"""
import sys

import numpy as np
import pysu2
from mpi4py import MPI

cfg, mode, nsteps = sys.argv[1], sys.argv[2], int(sys.argv[3])
amp = float(sys.argv[4]) if len(sys.argv) > 4 else 5e-5
assert mode in ("single", "reload-same", "reload-iqn") or mode.startswith("file:")
# file:<npz>[:nox|:noy] : apply a recorded fluid-interface displacement field (X, D) at every step
FIELD = None
if mode.startswith("file:"):
    parts = mode.split(":"); FIELD = dict(np.load(parts[1])); FLAG = parts[2] if len(parts) > 2 else ""
    if FLAG == "nox": FIELD["D"][:, 0] = 0.0
    if FLAG == "noy": FIELD["D"][:, 1] = 0.0
comm = MPI.COMM_WORLD
drv = pysu2.CSinglezoneDriver(cfg, 1, comm)

markers = drv.GetAllBoundaryMarkers()
mid = markers.get("interface")
verts, X = [], []
if mid is not None:
    for iv in range(drv.GetNumberVertices(mid)):
        if not drv.IsAHaloNode(mid, iv):
            verts.append(iv)
            X.append(drv.GetInitialMeshCoord(mid, iv)[:2])
x0, x1 = 0.2489898, 0.6


def set_disp(a):
    if FIELD is not None:
        for iv, (x, y) in zip(verts, X):
            k = int(np.argmin(np.hypot(FIELD["X"][:, 0] - x, FIELD["X"][:, 1] - y)))
            assert np.hypot(FIELD["X"][k, 0] - x, FIELD["X"][k, 1] - y) < 1e-5  # vtu export precision; nodes >= 3e-3 apart
            drv.SetMeshDisplacement(mid, iv, float(FIELD["D"][k, 0]), float(FIELD["D"][k, 1]), 0.0)
        return
    for iv, (x, y) in zip(verts, X):
        s = max(0.0, (x - x0) / (x1 - x0))
        drv.SetMeshDisplacement(mid, iv, 0.0, a * s * s, 0.0)


def solve(ti):
    drv.Preprocess(ti)
    drv.Run()
    drv.Postprocess()
    drv.Update()
    drv.Monitor(ti)


ti = drv.GetTime_Iter()
for n in range(1, nsteps + 1):
    a = amp * n
    if mode == "single" or FIELD is not None:
        set_disp(a)
        solve(ti)
    else:
        drv.SaveOldState()
        set_disp(a if mode == "reload-same" else 0.5 * a)
        solve(ti)
        drv.ReloadOldState()
        set_disp(a)
        solve(ti)
    drv.Output(ti)
    ti += 1
drv.Postprocessing()
if comm.Get_rank() == 0:
    print(f"done: mode {mode}, {nsteps} steps, amp {amp}")
