# Q1 - 2D left atrium of Qureshi et al. (CinC 2021)

## Calibration of the LAA wall velocity (CW)

| rhythm | W (mm/s) | Vmean LAA (cm/s) |
|---|---|---|
| SR | 0 | 0.320 |
| SR | **0.00** (chosen) | target 0.12 |
| AF | 0 | 0.163 |
| AF | **0.00** (chosen) | target 0.022 |

## Kinetics: qureshi

| shape | rhythm | seed | Vmean (cm/s) [Q] | Fn peak (uM) | Fn_res 60 s (uM) [Q] | t_TF (s) [Q] | t_washout (s) [Q] |
|---|---|---|---|---|---|---|---|
| CW | SR | P2 | 0.285 [0.12] | 4.50 | 4.50 [0.1] | 2.6 [-] | nan [28] |
| CW | AF | P2 | 0.188 [0.022] | 4.50 | 4.50 [3.7] | 2.6 [17] | nan [-] |

[Q] = Qureshi et al. 2021, Table 1 (seed P2 only).

Reference: Qureshi A, et al. Computing in Cardiology 48 (2021). doi:10.22489/CinC.2021.241
