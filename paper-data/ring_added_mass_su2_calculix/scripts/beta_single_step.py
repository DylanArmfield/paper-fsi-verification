#!/usr/bin/env python3
"""BETA_FACTOR check: one SU2 time step of the moderate-level, F1 fluid with a prescribed ring motion.

The interface moves by one time step of the initial n = 2 ovalling velocity (u_r = V0 dt cos 2 theta,
u_theta = -0.425 V0 dt sin 2 theta, the ring's outer-surface velocity), starting from the potential-flow
initial condition; prints the n = 2 modal interface force F2 = sum_i f_r,i cos 2 theta_i (N per m depth)
and the uniform part F0 = sum_i f_r,i.

Run in a folder holding fluid.su2 and restart_flow_00000/00001.csv (= restart_ic.csv) of
`setup_case.py moderate 1`, with ../beta/single_step.cfg (BETA_FACTOR, INNER_ITER as in beta/beta_single_step.csv):
    mpirun -np 2 python3 beta_single_step.py single_step.cfg > b.log
"""
import sys, math, pysu2
from mpi4py import MPI
d = pysu2.CSinglezoneDriver(sys.argv[1], 1, MPI.COMM_WORLD)
mid = d.GetAllBoundaryMarkers()["interface"]
V0, dt = 1.5e-3, d.GetUnsteady_TimeStep()
for iv in range(d.GetNumberVertices(mid)):
    x, y = d.GetInitialMeshCoord(mid, iv)[:2]; t = math.atan2(y, x)
    ur = V0*dt*math.cos(2*t); ut = V0*dt*math.sin(2*t)*(-0.425)
    d.SetMeshDisplacement(mid, iv, ur*math.cos(t)-ut*math.sin(t), ur*math.sin(t)+ut*math.cos(t), 0.0)
d.Preprocess(2); d.Run(); d.Postprocess(); d.Update(); d.Monitor(2)
F2 = F0 = 0.0
for iv in range(d.GetNumberVertices(mid)):
    x, y = d.GetInitialMeshCoord(mid, iv)[:2]; t = math.atan2(y, x)
    fx, fy = d.GetFlowLoad(mid, iv)[:2]
    fr = fx*math.cos(t) + fy*math.sin(t)
    F2 += fr*math.cos(2*t); F0 += fr
print("MODAL F2 %.10e  F0(uniform) %.10e" % (F2, F0))
