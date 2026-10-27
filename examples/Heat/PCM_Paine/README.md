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

Mesh: h_x = 5 mm, h_y = min(e/20, 1.25 mm). `REFINE=2` halves both for a
mesh-independence check. Environment variables of `run_pcm_sweep.sh` (ANGLE,
FACING, NDAYS, DT, UCEIL, ACH, ...) are listed at its top.

## Known issue: Newton in thick layers

Column tests (0.1 m x 50 mm, production h, dt 30 s): with the default 1 K
mushy range Newton diverges at the onset of melting (~9.7 h; dt 15 s fails
earlier, so it is not a time-step issue). With 2 K and `-snes_max_it 100` it
gets through melting (up to 63 iterations per step) but fails again during
solidification under the cooling top (~17 h, unstable layer). 10 mm passes the
onset of melting with the default settings (<= 22 iterations).
`run_pcm_sweep.sh` therefore restarts a failed case with 2 K +
`-snes_max_it 100`, then also with dt/2; `summary.md` records dt and the mushy
range actually used per case. A failed case keeps its logs as
`run_failed_dt*_mushy*.log`.
