# E4 - Patch-size threshold, no flow (Kastrup 2006)

| k_a | patch (um) | burst at (min) | gel at (min) | max T (nM) | front speed (um/min) | measured |
|---|---|---|---|---|---|---|
| 0.1 | 25 | 2.9 | 0.7 | 671.4 | 14.1 | no clot |
| 0.1 | 50 | 2.1 | 0.5 | 673.3 | 12.4 | no clot |
| 0.1 | 100 | 1.5 | 0.4 | 676.0 | 12.8 | clots < 3 min |
| 0.1 | 200 | 1.2 | 0.4 | 680.2 | 12.3 | clots < 3 min |
| 0.1 | 400 | 1.2 | 0.4 | 683.8 | 11.6 | clots < 3 min |
| 0.035 | 25 | nan | 0.7 | 14.4 | nan | no clot |
| 0.035 | 50 | 11.3 | 0.5 | 148.3 | 26.6 | no clot |
| 0.035 | 100 | 5.2 | 0.4 | 178.0 | 29.8 | clots < 3 min |
| 0.035 | 200 | 3.1 | 0.4 | 208.0 | 31.4 | clots < 3 min |
| 0.035 | 400 | 2.5 | 0.4 | 231.3 | 32.7 | clots < 3 min |

- k_a 0.1: threshold all clot (measured 50-100 um) -> FALLA
- k_a 0.035: threshold 25-400 um (measured 50-100 um) -> FALLA
- j_X = 4.0e-13 mol m^-2 s^-1 from Gamma_TF = 8e-13 mol/m^2 and kappa_TF = 0.5 1/s (E1), not fitted.

References:
- Kastrup CJ, Runyon MK, Shen F, Ismagilov RF. Modular chemical mechanism predicts spatiotemporal dynamics of initiation in the complex network of hemostasis. PNAS 103 (2006) 15747-15752
- Kastrup CJ, Shen F, Runyon MK, Ismagilov RF. Characterization of the threshold response of initiation of blood clotting to stimulus patch size. Biophys J 93 (2007) 2969-2977
