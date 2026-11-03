# PCM roof sandwich (Al / RT42 / Al): thickness sweep

`NaturalConvectionPCM` (enthalpy-porosity + natural convection, Newton + OSGS)
on a 2 m strip along the slope of a north-facing roof (25 deg by default):
RT42 paraffin between two 10 mm aluminium plates (`-pcm_al`, cells tagged
pcm / aluminium), ventilated attic underneath. Outside: solar gain, convection
h_a = 8.91 + 2 w and long-wave radiation to the sky; inside: convection only
(h_in = 2 W/m2K) to the attic air.

| File | Purpose |
|---|---|
| `paine_weather.py` | Clear January day, `--site paine` (-33.81, 370 m; GHI 8.2 kWh/m2/day, T_a 13-31 C, wind 1.5-4 m/s) or `--site calama` (-22.50, 2326 m; 8.3, 6-26 C, 2-9 m/s): irradiance on the roof plane (Hottel + isotropic sky) |
| `run_pcm_sweep.sh` | Runs e = 5, 10, 15, 20 mm of PCM (2 periodic days each; `SITE`, `AL`, `THICKNESSES`, ... at its top) |
| `postprocess_pcm.sh` | Heat-flux summary (`summary.md/csv`), plots and temperature contours |

```sh
cmake --build build -j --target NaturalConvectionPCM
mkdir -p ~/runs && cd ~/runs
NP=4 NJOBS=2 bash ~/rodin/examples/Heat/PCM_Paine/run_pcm_sweep.sh
bash ~/rodin/examples/Heat/PCM_Paine/postprocess_pcm.sh pcm_runs --hours 9 13 17 21
```

Bottom boundary (attic, "entretecho"): quasi-steady attic air node between the
PCM (h_in = 6 W/m2K, ISO 6946 downward), the ceiling to the room at 24 C
(U_c = 0.38 W/m2K, OGUC 4.1.10) and outdoor air (5 ACH, 0.9 m of attic per m2
of roof). Eliminating the node gives a Robin condition h_b (T - T_b(t)); see the
header of `NaturalConvectionPCM.cpp`. `-pcm_attic 0` restores the room-facing
condition h_in (T - T_in).

Mesh: h_x = 10 mm; 4 cells across each plate; across the PCM 0.25 mm at its
top face growing by 1.08 downwards (13 cells at 5 mm, 26 at 20 mm). `REFINE=2`
doubles the counts for a mesh-independence check. Environment variables of `run_pcm_sweep.sh` (ANGLE,
FACING, NDAYS, DT, UCEIL, ACH, ...) are listed at its top.

## Solver changes and checks (Nov 2026)

- **Aluminium sandwich** (`-pcm_al`, `-pcm_al_n`): cells tagged pcm /
  aluminium; plates conduct (rho c, k of aluminium), lumped enthalpy split by
  region. Flow only in the PCM: momentum and continuity integrated over it, the
  interfaces are no-slip walls (Dirichlet on interior facets tagged
  `interface`), OSGS projections restricted to the PCM. Checks: energy balance
  of the sandwich closes to 2e-4 over one hour; hydrostatic balance of the
  along-slope buoyancy exact (dp/dx 23.674 vs 23.678 Pa/m), zero net flow
  through any section, antisymmetric shear profile. (A Darcy wall plus flow
  equations in the plates gave spurious pressure modes and 5% unbalanced
  buoyancy, i.e. a spurious net flow.)
- **Exact Newton Jacobian**: A(T), H'(T), kappa(T) evaluated at the quadrature
  points of T_h (the P1 interpolants of nodal values had an inexact derivative
  where A' changes by orders of magnitude between nodes: linear convergence or
  2-cycles). Plus H''(T) dT u.grad T. Carman-Kozeny implicit in the residual,
  A'(T) dT u dropped from the Jacobian (`-pcm_darcy 2`, default; 0 exact,
  1 explicit A(T^n)). The line-search post-check halves updates that reverse
  near the mushy interval. Without plates and with the old settings the
  results match the previous version step by step.
- **Rodin SNES cache with PETSc >= 3.20** keyed on (object id, state): before,
  equal states of different vectors gave residuals at a stale state.
- **Sweep settings**: 2 K mushy range (RT42 melts over 38-43 C; with the
  solver default of 1 K, 15-20 mm lose Newton once the melt flows at
  ~3e-4 m/s along the slope) and up to 100 Newton iterations. Column tests
  (0.1 m, Calama): 5 and 10 mm full day; 15 and 20 mm from 09:00 to 17:00
  through melting, peak and onset of solidification (<= 13 iterations per
  step). A failed case restarts with 4 K; `summary.md` records the mushy range
  of each case. A smaller dt does not help (fails earlier, same state).
