#!/usr/bin/env python3
"""Static check of the Kratos plate model, the same diagnostic as solids4foam's solid_discretisation.py:
the plate of beam_case.py alone, clamped at y = 0, symmetric at z = 0, a uniform traction of 100 Pa
in +x on the upstream face x = 0.4 and the other faces free. u_x at A against solids4foam's static
sequence (paper-data/beam_richter/solid_static.csv, family m2: 8f cells through the thickness, 8f in
y and z), which converges to about 1.41e-4 m.

Usage: python scripts/solid_static_check.py <work dir> <factor> [...]
         [--element TotalLagrangianElement3D8N]
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beam_case as bc  # noqa: E402

TRACTION = 100.0
S4F_M2 = {1: 1.23944e-4, 2: 1.35100e-4, 4: 1.39520e-4, 8: 1.40874e-4}


def write_case(work, f, element):
    os.makedirs(work)
    info = bc.write_solid(os.path.join(work, bc.SOLID_MDPA + ".mdpa"), f, 2, 0, 0, element=element,
                          condition="SurfaceLoadCondition3D3N", interface_faces=("x0",))
    bc.write_materials(work)
    params = {
        "problem_data": {"problem_name": "plate", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": 1.0},
        "solver_settings": {"solver_type": "Static", "model_part_name": "Structure", "domain_size": 3,
                            "echo_level": 0, "analysis_type": "non_linear",
                            "model_import_settings": {"input_type": "mdpa", "input_filename": bc.SOLID_MDPA},
                            "material_import_settings": {"materials_filename": "StructuralMaterials.json"},
                            "time_stepping": {"time_step": 1.0}, "rotation_dofs": False,
                            "convergence_criterion": "residual_criterion",
                            "residual_relative_tolerance": 1e-10, "residual_absolute_tolerance": 1e-12,
                            "max_iteration": 20},
        "processes": {
            "constraints_process_list": [
                bc.vector_process("Structure.Clamped", "DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                bc.vector_process("Structure.Symmetry", "DISPLACEMENT", [False, False, True], [None, None, 0.0])],
            "loads_process_list": [{
                "python_module": "assign_vector_by_direction_to_condition_process",
                "kratos_module": "KratosMultiphysics", "process_name": "AssignVectorByDirectionToConditionProcess",
                "Parameters": {"model_part_name": "Structure.StructureInterface", "variable_name": "SURFACE_LOAD",
                               "modulus": TRACTION, "direction": [1.0, 0.0, 0.0], "interval": [0.0, "End"]}}]},
        "output_processes": {}}
    with open(os.path.join(work, "ProjectParameters.json"), "w") as fo:
        json.dump(params, fo, indent=2)
    return info


def run(work):
    import KratosMultiphysics as KM
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
    cwd = os.getcwd()
    os.chdir(work)
    with open("ProjectParameters.json") as f:
        params = KM.Parameters(f.read())
    model = KM.Model()
    StructuralMechanicsAnalysis(model, params).Run()
    mp = model.GetModelPart("Structure")
    a = bc.POINT_A
    node = min(mp.Nodes, key=lambda n: (n.X0 - a[0])**2 + (n.Y0 - a[1])**2 + (n.Z0 - a[2])**2)
    assert math.dist((node.X0, node.Y0, node.Z0), a) < 1e-10
    out = [node.GetSolutionStepValue(v) for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)]
    load = sum(c.GetGeometry().Area() for c in model.GetModelPart("Structure.StructureInterface").Conditions) * TRACTION
    os.chdir(cwd)
    return out, load


def main():
    argv = list(sys.argv[1:])
    element = "TotalLagrangianElement3D8N"
    if "--element" in argv:
        i = argv.index("--element")
        element = argv[i + 1]
        del argv[i:i + 2]
    work_root, factors = argv[0], [int(f) for f in argv[1:]]
    rows = []
    for f in factors:
        work = os.path.join(work_root, f"plate_{element}_f{f}")
        info = write_case(work, f, element)
        u, load = run(work)
        ref = S4F_M2.get(f)
        rows.append({"factor": f, "element": element, "cells_xyz": info["solid_cells_xyz"], "load_N": load,
                     "ux_A": u[0], "uy_A": u[1], "uz_A": u[2], "s4f_m2_ux_A": ref})
        print(f"factor {f} {info['solid_cells_xyz']} {element}: load {load:.6f} N, u_x(A) {u[0]:.6e} m"
              + (f" (solids4foam m2 {ref:.6e}, {100 * (u[0] / ref - 1):+.2f}%;" if ref else " (")
              + f" vs 1.41e-4: {100 * (u[0] / 1.41e-4 - 1):+.2f}%), u_y {u[1]:.4e}, u_z {u[2]:.4e}")
    with open(os.path.join(work_root, f"solid_static_{element}.json"), "w") as fo:
        json.dump(rows, fo, indent=2)


if __name__ == "__main__":
    main()
