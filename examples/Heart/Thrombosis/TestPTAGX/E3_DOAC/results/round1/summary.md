# E3 - Thrombodynamics under DOACs, PTAGX kinetics (k_a = 0.035)

| case | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) | Tlag/ctrl | Vi/ctrl | CS30/ctrl |
|---|---|---|---|---|---|---|---|
| control | 0.75 | 39.6 | 28.8 | 1085 | 1.00 | 1.00 | 1.00 |
| Xa_both_0.97 | 0.75 | 38.6 | 20.8 | 920 | 1.00 | 0.97 | 0.85 |
| Xa_both_0.94 | 0.83 | 36.7 | 14.8 | 790 | 1.11 | 0.93 | 0.73 |
| Xa_both_0.9 | 0.83 | 35.6 | 10.3 | 675 | 1.11 | 0.90 | 0.62 |
| Xa_both_0.85 | 0.83 | 34.2 | 7.3 | 570 | 1.11 | 0.86 | 0.53 |
| Xa_both_0.8 | 0.83 | 33.2 | 5.2 | 500 | 1.11 | 0.84 | 0.46 |
| Xa_kappa_0.5 | 1.08 | 30.4 | 32.8 | 930 | 1.44 | 0.77 | 0.86 |
| Xa_kappa_0.25 | 1.50 | 24.4 | 3.7 | 360 | 2.00 | 0.62 | 0.33 |
| Xa_kappa_0.1 | 2.33 | 18.7 | 3.7 | 290 | 3.11 | 0.47 | 0.27 |
| IIa_1.05 | 0.75 | 38.0 | 13.2 | 795 | 1.00 | 0.96 | 0.73 |
| IIa_1.1 | 0.75 | 36.0 | 6.2 | 645 | 1.00 | 0.91 | 0.59 |
| IIa_1.2 | 0.75 | 33.3 | 2.8 | 490 | 1.00 | 0.84 | 0.45 |
| IIa_1.35 | 0.83 | 29.7 | 3.1 | 425 | 1.11 | 0.75 | 0.39 |

Peak/trough ratios: measured vs model at the dose that matches the Tlag ratio.

| drug | family | measured Tlag/Vi/CS30 | model dose | model Tlag/Vi/CS30 | Vi, CS30 within 0.05 |
|---|---|---|---|---|---|
| rivaroxaban | Xa_both | 1.67/0.86/0.73 | Tlag ratio not reached (model 1.00-1.11) | - | NO |
| rivaroxaban | Xa_kappa | 1.67/0.86/0.73 | between Xa_kappa_0.5 and Xa_kappa_0.25 | 1.67/0.71/0.65 | NO |
| apixaban | Xa_both | 1.18/0.95/0.82 | Tlag ratio not reached (model 1.00-1.11) | - | NO |
| apixaban | Xa_kappa | 1.18/0.95/0.82 | between control and Xa_kappa_0.5 | 1.18/0.91/0.94 | NO |
| dabigatran | IIa | 1.30/0.87/0.84 | Tlag ratio not reached (model 1.00-1.11) | - | NO |

References:
- ADAGE substudy: Fibrinography and thrombography (Thrombodynamics-4D) in atrial fibrillation: assessment of direct oral anticoagulants in geriatric patients aged 80 years and older (PMC12351382), Tables 1-3
- Gebekova ZA, et al. Ration Pharmacother Cardiol 18 (2022) 544-552
