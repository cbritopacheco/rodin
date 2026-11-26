# E3 - Thrombodynamics under DOACs, PTAGX kinetics (MODEL)

| case | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) | Tlag/ctrl | Vi/ctrl | CS30/ctrl |
|---|---|---|---|---|---|---|---|
| control | 1.00 | 55.4 | 26.8 | 1160 | 1.00 | 1.00 | 1.00 |
| Xa_both_0.97 | 1.00 | 53.1 | 17.0 | 975 | 1.00 | 0.96 | 0.84 |
| Xa_both_0.94 | 1.00 | 50.7 | 10.3 | 850 | 1.00 | 0.92 | 0.73 |
| Xa_both_0.9 | 1.08 | 48.0 | 6.0 | 745 | 1.08 | 0.87 | 0.64 |
| Xa_both_0.85 | 1.08 | 45.2 | 3.7 | 655 | 1.08 | 0.82 | 0.56 |
| Xa_both_0.8 | 1.08 | 42.7 | 3.1 | 595 | 1.08 | 0.77 | 0.51 |
| Xa_kappa_0.9 | 1.08 | 53.6 | 27.2 | 1155 | 1.08 | 0.97 | 1.00 |
| Xa_kappa_0.75 | 1.17 | 50.0 | 27.9 | 1140 | 1.17 | 0.90 | 0.98 |
| Xa_kappa_0.6 | 1.33 | 45.3 | 28.9 | 1115 | 1.33 | 0.82 | 0.96 |
| Xa_kappa_0.5 | 1.42 | 40.9 | 30.0 | 1085 | 1.42 | 0.74 | 0.94 |
| Xa_kappa_0.25 | 2.17 | 22.9 | 34.8 | 800 | 2.17 | 0.41 | 0.69 |
| IIa_1.05 | 1.00 | 52.1 | 9.1 | 850 | 1.00 | 0.94 | 0.73 |
| IIa_1.1 | 1.00 | 49.7 | 2.9 | 705 | 1.00 | 0.90 | 0.61 |
| IIa_1.2 | 1.00 | 45.3 | 0.8 | 560 | 1.00 | 0.82 | 0.48 |
| IIa_1.35 | 1.00 | 39.5 | 0.7 | 455 | 1.00 | 0.71 | 0.39 |

Peak/trough ratios: measured vs model at the dose that matches the Tlag ratio.

| drug | family | measured Tlag/Vi/CS30 | model dose | model Tlag/Vi/CS30 | Vi, CS30 within 0.05 |
|---|---|---|---|---|---|
| rivaroxaban | Xa_both | 1.67/0.86/0.73 | Tlag ratio not reached (model 1.00-1.08) | - | NO |
| rivaroxaban | Xa_kappa | 1.67/0.86/0.73 | between Xa_kappa_0.5 and Xa_kappa_0.25 | 1.67/0.63/0.85 | NO |
| apixaban | Xa_both | 1.18/0.95/0.82 | Tlag ratio not reached (model 1.00-1.08) | - | NO |
| apixaban | Xa_kappa | 1.18/0.95/0.82 | between Xa_kappa_0.75 and Xa_kappa_0.6 | 1.18/0.90/0.98 | NO |
| dabigatran | IIa | 1.30/0.87/0.84 | Tlag ratio not reached (model 1.00-1.00) | - | NO |

References:
- ADAGE substudy: Fibrinography and thrombography (Thrombodynamics-4D) in atrial fibrillation: assessment of direct oral anticoagulants in geriatric patients aged 80 years and older (PMC12351382), Tables 1-3
- Gebekova ZA, et al. Ration Pharmacother Cardiol 18 (2022) 544-552
