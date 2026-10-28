# PCM roof layer under Paine weather: thickness sweep

`NaturalConvectionPCM` (enthalpy-porosity + natural convection, Newton + OSGS)
on a 2 m strip along the slope of a north-facing roof (25 deg by default),
RT42 paraffin, ventilated attic underneath.

| File | Purpose |
|---|---|
| `paine_weather.py` | Clear January day in Paine (-33.81, 370 m): irradiance on the roof plane (Hottel + isotropic sky, GHI scaled to 8.2 kWh/m2/day), T_a 13-31 C, valley wind 1.5-4 m/s |
| `run_pcm_sweep.sh` | Runs e = 5, 10, ..., 50 mm (3 periodic days each) |
| `postprocess_pcm.sh` | Heat-flux summary (`summary.md/csv`), plots and temperature contours |

```sh
cmake --build build -j --target NaturalConvectionPCM
mkdir -p ~/runs && cd ~/runs
NP=4 NJOBS=2 bash ~/rodin/examples/Heat/PCM_Paine/run_pcm_sweep.sh
bash ~/rodin/examples/Heat/PCM_Paine/postprocess_pcm.sh pcm_paine_runs --hours 9 13 17 21
```

Bottom boundary (attic, "entretecho"): quasi-steady attic air node between the
PCM (h_in = 6 W/m2K, ISO 6946 downward), the ceiling to the room at 24 C
(U_c = 0.38 W/m2K, OGUC 4.1.10) and outdoor air (5 ACH, 0.9 m of attic per m2
of roof). Eliminating the node gives a Robin condition h_b (T - T_b(t)); see the
header of `NaturalConvectionPCM.cpp`. `-pcm_attic 0` restores the room-facing
condition h_in (T - T_in).

Mesh: h_x = 5 mm; across the layer 0.25 mm at the top face growing by 1.04
downwards (20 cells at 5 mm, 57 at 50 mm). `REFINE=2` doubles both counts for
a mesh-independence check. Environment variables of `run_pcm_sweep.sh` (ANGLE,
FACING, NDAYS, DT, UCEIL, ACH, ...) are listed at its top.

## Newton robustness (findings, 28 Oct 2026)

- **Mesh at the top face.** The first liquid layer must span several cells.
  Column tests (0.1 m x 50 mm, 25 deg, dt 30 s): with a uniform h_y of 1.25 mm
  Newton diverges ~1.5 h after the onset of melting (iterations 5 -> 30);
  0.5 mm also climbs; 0.25 mm keeps ~9 iterations past that point. The sweep
  grades the mesh from 0.25 mm at the top (`-pcm_ygrade`), which delays but
  does not remove the failure at 50 mm (diverges at Y = 0.065, liquid layer
  ~3 mm, when along-slope flow develops).
- **Open: thick layers on a slope.** The along-slope flow in the melt grows
  with the layer (max |u| 4e-5 m/s at 5 mm, the run that worked; 1e-4 to
  2e-4 m/s at 10-50 mm) and Newton degrades with it: 10 mm (graded) diverged
  at Y = 0.68, 50 mm at Y = 0.065. Not cured by a smaller dt (15 s fails
  earlier), the full H'' Jacobian term, or backtracking line search; a 2 K
  mushy range with `-snes_max_it 100` carried 50 mm through melting but it
  failed at ~17 h. Needs algorithmic work (e.g. step-level dt retries from
  t^n, or Picard iterations before Newton).
- **Rodin SNES cache with PETSc >= 3.20.** `Rodin::Solver::SNES` skipped the
  state update / reassembly when `VecGetState(x)` matched the last one, but
  states are per object and SNES evaluates the residual at the line-search
  vector W and the Jacobian at X: equal states with different contents gave
  residuals at a stale state (reproduced, |x - x_state| = 1e-3). The cache key
  is now (object id, state); with it the Newton history is identical step by
  step to the uncached PETSc 3.19 path.
- If a case still fails, `run_pcm_sweep.sh` restarts it with a 2 K mushy range
  (RT42 melts over 38-43 C) and `-snes_max_it 100`, then also dt/2;
  `summary.md` records dt and the mushy range used, and the failed logs are
  kept as `run_failed_dt*_mushy*.log`.
