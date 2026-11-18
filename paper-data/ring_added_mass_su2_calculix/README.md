# ringAddedMass with preCICE: SU2 fluid + CalculiX solid (November 2026)

This is an independent-code solution of the ring added-mass problem (issue #36). The fluid code (SU2, a vertex-based finite volume solver), the solid code (CalculiX, finite elements) and the coupling library (preCICE, IQN-ILS) are all independent of solids4foam. The route is partitioned, so it complements FEBio's monolithic solution in `../febio_ring_added_mass`.

The problem and the exact frequencies are those of the solids4foam `ringAddedMass` verification: its README and reference JSON are copied in `../febio_ring_added_mass/reference/s4f_ringAddedMass/verification/`. The scripts here read the JSON from there. The platform and the SU2 build are those of `../fsi1_su2_calculix` (README, Platform), and the environment script is `../fsi1_su2_calculix/scripts/env-su2-precice.sh`.

## Set-up

- **Geometry and material:** ring 0.95–1.05 m inside a rigid circular wall at 1.5 m, plane strain. E = 1 MPa, ν = 0.3, ρs = 1000 kg/m³. ρf = 14 / 140 / 1400 kg/m³ for the weak / moderate / strong levels. The full ring is modelled, not a quadrant.
- **Meshes** (`scripts/make_ring.py`): for refinement factor F = 1, 2, 4:
  - fluid: 12F × 48F quadrilaterals, uniform in r and θ;
  - solid: 4F × 48F CalculiX C3D8I elements, one layer of 1 m depth with u_z = 0;
  - the interface nodes coincide (nearest-neighbour mapping).
- **Initial condition:** the n = 2 ovalling velocity with V0 = 1.5e-3 m/s and zero displacement in the ring (as the FEBio evidence), with the matching potential flow in the fluid and p = 0. The fluid starts from it through an ASCII SU2 restart, with both BDF2 levels set to it.
- **Fluid** (`case/fluid-su2/coupled.cfg`):
  - SU2 7.5.1 incompressible Euler (`INC_EULER`), slip on the interface and the outer wall: the reference is the inviscid solution;
  - FDS with MUSCL reconstruction and no limiter, weighted-least-squares gradients, `BETA_FACTOR` 1e6 (see the β section);
  - BDF2 dual time stepping with Δt = the coupling window. Pseudo-time with implicit Euler, CFL 1000, FGMRES + ILU (1e-8), 300 inner iterations per solve;
  - mesh motion by linear elasticity with inverse-volume stiffness.
- **Solid** (`case/solid-calculix/ring.inp`): CalculiX 2.20 (`ccx_preCICE`), linear elastic and linear geometry (no `NLGEOM`), `*DYNAMIC` with HHT α = −0.05. The **dry** runs use the same binary and input without coupling.
- **Coupling** (`case/precice-config.base.xml`):
  - preCICE 3.5.0 parallel-implicit, with a relative convergence limit of 1e-6 on force and displacement and at most 50 iterations;
  - IQN-ILS with residual-sum preconditioner, QR2 filter 1e-3, initial relaxation 0.1, 60 used iterations;
  - **10 reused time windows at the weak and moderate levels and 0 at the strong level** (see Coupling).
- **Time stepping:** 100 steps per period of the exact frequency of each level (as the solids4foam driver), 4 periods.
- **Evaluation** (`scripts/evaluate_ring.py`, as the solids4foam verification driver):
  - The signal is q(t) = (u_r(0°) − u_r(90°))/2 on the outer surface. It is taken from the CalculiX node prints at (1.05, 0) and (0, 1.05).
  - The frequency comes from the zero crossings, interpolated linearly: the half-period is the least-squares slope of crossing time against crossing index.
  - The amplitude loss per period is measured from the peaks, with parabolic interpolation.
  - The wet/dry ratio divides the wet frequency by that of a dry run on the same mesh and step.

## Mesh study

Wet/dry ratio error against the exact continuum ratio, at 100 steps per period (`ring_mesh_study.csv`):

| Level | F1 | F2 | F4 | Observed order | Richardson ratio error | Frequency error at F4 | Loss per period at F4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| weak (m* = 0.10) | +0.104 % | +0.030 % | **+0.010 %** | 1.92 | +0.003 % | −0.049 % | 0.007 % |
| moderate (m* = 1.0) | +0.557 % | +0.150 % | **+0.042 %** | 1.92 | +0.004 % | −0.017 % | 0.029 % |
| strong (m* = 10) | +0.986 % | +0.247 % | **+0.052 %** | 1.91 | −0.019 % | −0.008 % | 0.007 % |

- **Ratio convergence:** the ratio error converges at second order at every level and extrapolates to within 0.02 % of the exact ratio.
- **The dry ring:** its frequency error is +0.676 / +0.087 / −0.059 % at F1 / F2 / F4. The C3D8I ring is therefore not converged in space on the coarse meshes, as with solids4foam's standard solid, and the ratio removes most of this common error.
- **Other codes at F4 or their finest mesh:**
  - solids4foam (high-order solid, IQN-ILS): +0.001 / +0.003 / +0.006 % at F4;
  - FEBio: +0.004 / +0.008 / +0.041 % at its finest mesh, with a floor that is nearly independent of the mesh.

## Time-step study

Moderate level and the dry ring, F2 (`ring_time_study.csv`):

| | 25 | 50 | 100 steps/period | Observed order |
|---|---:|---:|---:|---:|
| dry frequency error | −0.469 % | −0.026 % | +0.087 % | 1.98 |
| moderate frequency error | −0.975 % | −0.015 % | +0.237 % | 1.93 |
| moderate ratio error | −0.508 % | +0.011 % | +0.150 % | |
| moderate loss per period | 1.39 % | 0.36 % | 0.24 % | |

- **Time error:** it is second order. At 100 steps per period it is about −0.08 % in the frequency and about −0.05 % in the ratio, so the F2 ratio error of +0.15 % is mostly spatial.
- **The near-zero ratio error at 50 steps** is a cancellation of the spatial and time errors, which have opposite signs.

## Acceptance criteria

The solids4foam verification criteria, with the tolerances from the reference JSON (`ring_acceptance.csv`):

- **Met:**
  - ratio error at the finest mesh within 0.1 % at all three levels, decreasing under refinement;
  - observed mesh order of the ratio at least 1.5 (1.91–1.92);
  - ratio error against thin-ring theory at the finest mesh within 0.25 % (+0.02 / +0.11 / +0.17 % for weak / moderate / strong);
  - frequency error at the finest mesh within 0.3 % (dry and all three levels);
  - amplitude loss at F2 below 2 % per period (0.05 / 0.24 / 0.52 %), with no run growing;
  - observed time order at least 1.5 (1.98 dry, 1.93 moderate).
- **Not met, because the C3D8I solid is not converged in space:**
  - the time-extrapolated frequency within 0.1 %: +0.12 % (dry) and +0.32 % (moderate) at F2. The extrapolation removes the time error only, so these values contain the F2 spatial error of the ring.
  - the frequency change between meshes within 0.1 %: from F2 to F4 the frequency changes by −0.15 % (dry) and −0.17 / −0.25 / −0.34 % (weak / moderate / strong).
  - Both criteria assume a solid that is converged in space on every mesh (the solids4foam high-order solid). solids4foam's standard solid is checked on the mesh order and the ratio instead, which are met here.
- **Not tested:** the Robin–Neumann against IQN-ILS agreement, because only IQN-ILS was run.

## SU2 `BETA_FACTOR`

SU2's incompressible solver uses artificial compressibility with β² = `BETA_FACTOR` · max|u|². Here max|u| is the 1.5 mm/s initial ring velocity, so the default β is tiny and the pseudo-time solve stalls.

- **Single-step start-up check** (`beta/`, `scripts/beta_single_step.py`, moderate level, F1). It was used to choose a β at which the SU2 solve converges; it is not evidence about the added mass.
  - the interface is moved by one time step of the initial ring velocity from the potential-flow state, and the n = 2 modal interface force is printed (`beta/beta_single_step.csv`);
  - at the default β (4.1) the pressure residual stalls (log10 rms[p] −2.6 / −3.0 / −4.0 after 50 / 300 / 2000 inner iterations), and the modal force is 1 / 5 / 7 % of the value at β 1e6, so the default β cannot be used;
  - from β 1e5 upwards, rms[p] reaches its floor of −6.96 (closed annulus: the pressure level floats). The single-step force still changes with β: 0.96 / 1.00 / 1.02 / 1.03 / 1.21 times the β 1e6 value at β 1e5 / 1e6 / 1e7 / 1e8 / 1e10. This change, about 4 % between β 1e5 and 1e8, does not appear in the coupled ratio, which moves by 0.02 % (below). The single-step force is therefore not a measure of the added mass.
- **Coupled runs** (moderate level, 100 steps per period). At F1 the β runs differ in the residual the inner solves reach, as well as in β, so F1 was repeated with 1000 inner iterations instead of 300, and β 1e5 / 1e8 were also run at F2:

  | Run | β 1e5 | β 1e6 | β 1e8 |
  |---|---:|---:|---:|
  | F1, 300 inner iterations: ratio error | +0.571 % | +0.557 % | +0.553 % |
  | F1, 300 inner iterations: loss per period | 3.47 % | 2.28 % | 10.02 % |
  | F1, 300: median final log10 rms[p], last coupling iteration | −6.71 | −7.48 | −8.58 |
  | F1, 1000 inner iterations: ratio error | +0.568 % | +0.557 % | +0.553 % |
  | F1, 1000 inner iterations: loss per period | 3.37 % | 2.26 % | 10.02 % |
  | F1, 1000: median final log10 rms[p], last coupling iteration | −7.51 | −7.99 | −9.10 |
  | F2, 300 inner iterations: ratio error | +0.149 % | +0.150 % | +0.162 % |
  | F2, 300 inner iterations: loss per period | 0.44 % | 0.24 % | 0.91 % |

  - **Inner iterations:** 1000 inner iterations lower the residual by 0.5–0.8 orders, but change the loss by at most 0.1 % per period and the ratio by at most 0.003 %. The β effect on the damping is therefore not an effect of the 300-iteration cap. All solves still stop at the cap.
  - **Mesh refinement:** from F1 to F2 the loss falls by a factor 8–11 at every β, and the spread across β falls from 7.7 to 0.67 % per period. At β 1e6 the loss falls again to 0.03 % at F4. The β dependence of the damping is a discretisation error that decreases under refinement, but at F2 it has not disappeared: β 1e8 still damps about four times as much as β 1e6.
  - **Ratio:** β moves the ratio by 0.02 % at F1 and 0.01 % at F2.
- **Conclusion:** as in FSI1 and CFD3 (`../fsi1_su2_calculix`), β is part of the discretisation. All the production runs use β 1e6.

## Coupling

- **Coupling iterations:** the mean per window at F4 is 7.2 (weak), 10.2 (moderate) and 22.9 (strong). No window reached the 50-iteration cap.
- **IQN-ILS reuse at the strong level:** with 10 reused windows (`strong-F1-reuse10`), the coupling iterations rose from about 22 to 34 per window and SU2 diverged in window 130. Without reuse the strong level runs at every mesh.
- **SU2 inner iterations:** every SU2 solve uses all its inner iterations (300 in the production runs, 1000 in the `-inner1000` runs), because the absolute stopping criterion (log10 rms −14) is below the pressure floor. The median final log10 rms[p] of the solves is −7.8 / −8.3 / −8.7 (weak), −7.0 / −7.4 / −7.9 (moderate) and −6.6 / −7.0 / −7.5 (strong) at F1 / F2 / F4. The coupling convergence is measured by preCICE on the interface data.

## Cost

Wall clock (`ring_su2_calculix.csv`), with cores given as SU2 ranks + 1 CalculiX:

| Level | F1 (2 + 1) | F2 (4 + 1) | F4 (8 + 1) |
|---|---:|---:|---:|
| weak | 9 min | 26 min | 1.9 h |
| moderate | 12 min | 30 min | 2.4 h |
| strong | 22 min | 1.2 h | 3.8 h |

Most of the cost is the 300 inner iterations per SU2 solve, multiplied by the coupling iterations. The dry runs take 6 to 127 s on one core.

## Files

- **Summary tables** (`scripts/ring_table.py` rebuilds them from `histories/` and `logs/` only):
  - `ring_su2_calculix.csv`: one row per run;
  - `ring_mesh_study.csv`, `ring_time_study.csv`, `ring_acceptance.csv`.
- **`histories/<run>.csv`:** per step, the time, u_r at 0° and 90° (as printed by CalculiX) and the coupling iterations of the window.
- **`histories/<run>_su2.csv`:** one row per SU2 solve, with the window, the coupling iteration, the inner iterations and the final log10 rms residuals.
- **`logs/<run>/`** for each run:
  - Slurm output, `CASE_INFO.txt`, `case.json`, exit status and wall time;
  - the preCICE config used and the preCICE iteration and convergence logs;
  - the SU2 `coupled.cfg`;
  - the first 120 and last 60 lines of the SU2 and CalculiX screen logs.
- **`case/`:** the templates filled in by `scripts/setup_case.py`.
- **`scripts/`:**
  - `make_ring.py`: meshes and initial condition;
  - `setup_case.py`: creates a case;
  - `run_ring.sh`: Slurm, dry or coupled, with a hang guard and a progress watchdog;
  - `evaluate_ring.py`: reads a case folder or a published history;
  - `ring_table.py`;
  - `beta_single_step.py`.
- **`beta/`:** the single-step cfg, the SU2 logs and `beta_single_step.csv`.

**Runs:**
- `dry-F{1,2,4}` and `dry-F2-spp{25,50}`;
- `{weak,moderate,strong}-F{1,2,4}`;
- `moderate-F2-spp{25,50}` and `moderate-F1-beta{1e5,1e8}`;
- `strong-F1-reuse10`: stopped at window 130, exit 1;
- β and inner-iteration runs: `moderate-F1-inner1000`, `moderate-F1-beta{1e5,1e8}-inner1000` and `moderate-F2-beta{1e5,1e8}`.

## Reproduction

From this folder:

```
W="--beta 1e6 --inner-iter 300 --inner-minval -14"
python3 scripts/setup_case.py dry 2 cases/dry-F2                    # dry: add --spp 25 / 50 for the time study
python3 scripts/setup_case.py moderate 2 cases/moderate-F2 $W       # weak / moderate; --spp 25 / 50; --beta 1e5 / 1e8; --inner-iter 1000 (after $W)
python3 scripts/setup_case.py strong 2 cases/strong-F2 $W --reuse 0
sbatch --ntasks=1 --export=ALL,CASE=cases/dry-F2,STAGE=dry scripts/run_ring.sh
sbatch --ntasks=5 --export=ALL,CASE=cases/moderate-F2,STAGE=coupled scripts/run_ring.sh   # F1 3, F2 5, F4 9 tasks
python3 scripts/evaluate_ring.py cases/moderate-F2 --dry cases/dry-F2
python3 scripts/evaluate_ring.py histories/moderate-F2.csv --dry histories/dry-F2.csv     # the same from the published history
python3 scripts/ring_table.py
```

**Checks made:**
- `setup_case.py` regenerates the inputs of all 24 runs byte for byte. At the strong level, IQN reuse was set to 0 by editing the generated config (now `--reuse 0`), so one XML comment and the `CASE_INFO.txt` wording differ.
- Evaluating the published histories gives the case-folder results to the 7 digits of the tables.

## Limitations

- **Iterative error:** the SU2 inner solves stop at the iteration cap, not at a residual criterion. At moderate F1, raising the cap from 300 to 1000 changed the ratio by at most 0.003 % and the loss by at most 0.1 % per period; the coupling tolerance was not varied.
- **β:** the frequency and ratio depend little on β (0.02 % at F1, 0.01 % at F2 between β 1e5 and 1e8). The amplitude loss depends on β; the dependence decreases under refinement but is still a factor of about four between β 1e6 and 1e8 at F2.
- **Coverage:**
  - the time-step study covers the moderate level only, at F2;
  - only IQN-ILS was run;
  - the amplitude dependence that FEBio found at the strong level (`../febio_ring_added_mass`, amplitude diagnostic) was not tested;
  - the finite-viscosity problem is not solved (inviscid fluid, slip walls).
