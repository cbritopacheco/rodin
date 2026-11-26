# E3 - Thrombodynamics under DOACs, PTAGX kinetics (MODEL)

| case | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) | Tlag/ctrl | Vi/ctrl | CS30/ctrl |
|---|---|---|---|---|---|---|---|
| control | 0.33 | 61.7 | 18.1 | 1165 | 1.00 | 1.00 | 1.00 |
| Xa_both_0.97 | 0.33 | 59.7 | 5.2 | 935 | 1.00 | 0.97 | 0.80 |
| Xa_both_0.94 | 0.33 | 58.1 | 1.6 | 825 | 1.00 | 0.94 | 0.71 |
| Xa_both_0.9 | 0.33 | 56.1 | 0.5 | 740 | 1.00 | 0.91 | 0.64 |
| Xa_both_0.85 | 0.33 | 53.0 | 0.4 | 675 | 1.00 | 0.86 | 0.58 |
| Xa_both_0.8 | 0.42 | 49.4 | 0.3 | 630 | 1.25 | 0.80 | 0.54 |
| Xa_kappa_0.5 | 0.42 | 61.7 | 18.2 | 1160 | 1.25 | 1.00 | 1.00 |
| Xa_kappa_0.25 | 0.58 | 61.1 | 18.4 | 1160 | 1.75 | 0.99 | 1.00 |
| Xa_kappa_0.1 | 0.83 | 61.3 | 18.8 | 1160 | 2.50 | 0.99 | 1.00 |
| IIa_1.05 | 0.33 | 58.1 | 1.3 | 820 | 1.00 | 0.94 | 0.70 |
| IIa_1.1 | 0.33 | 54.5 | 0.1 | 695 | 1.00 | 0.88 | 0.60 |
| IIa_1.2 | 0.33 | 46.7 | 0.1 | 575 | 1.00 | 0.76 | 0.49 |
| IIa_1.35 | 0.33 | 36.9 | 0.7 | 495 | 1.00 | 0.60 | 0.42 |

Peak/trough ratios: measured vs model at the dose that matches the Tlag ratio.

| drug | family | measured Tlag/Vi/CS30 | model dose | model Tlag/Vi/CS30 | Vi, CS30 within 0.05 |
|---|---|---|---|---|---|
| rivaroxaban | Xa_both | 1.67/0.86/0.73 | Tlag ratio not reached (model 1.00-1.25) | - | NO |
| rivaroxaban | Xa_kappa | 1.67/0.86/0.73 | between Xa_kappa_0.5 and Xa_kappa_0.25 | 1.67/0.99/1.00 | NO |
| apixaban | Xa_both | 1.18/0.95/0.82 | Tlag ratio not reached (model 1.00-1.25) | - | NO |
| apixaban | Xa_kappa | 1.18/0.95/0.82 | between control and Xa_kappa_0.5 | 1.18/1.00/1.00 | NO |
| dabigatran | IIa | 1.30/0.87/0.84 | Tlag ratio not reached (model 1.00-1.00) | - | NO |

References:
- ADAGE substudy: Fibrinography and thrombography (Thrombodynamics-4D) in atrial fibrillation: assessment of direct oral anticoagulants in geriatric patients aged 80 years and older (PMC12351382), Tables 1-3
- Gebekova ZA, et al. Ration Pharmacother Cardiol 18 (2022) 544-552
