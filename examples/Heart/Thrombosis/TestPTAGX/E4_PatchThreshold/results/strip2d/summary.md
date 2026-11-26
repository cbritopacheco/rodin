# E4 - Patch-size threshold, no flow (Kastrup 2006)

| set | patch (um) | burst at (min) | gel at (min) | max T (nM) | front speed (um/min) | measured |
|---|---|---|---|---|---|---|
| model | 25 | 0.5 | 0.4 | 1208.1 | 30.3 | no clot |
| model | 50 | 0.3 | 0.2 | 1256.2 | 30.2 | no clot |
| model | 100 | 0.2 | 0.2 | 1272.9 | 30.3 | clots < 3 min |
| model | 200 | 0.2 | 0.2 | 1274.8 | 30.3 | clots < 3 min |
| model | 400 | 0.2 | 0.2 | 1274.7 | 31.0 | clots < 3 min |
| legacy | 25 | 2.9 | 0.7 | 671.4 | 14.1 | no clot |
| legacy | 50 | 2.1 | 0.5 | 673.3 | 12.4 | no clot |
| legacy | 100 | 1.5 | 0.4 | 676.0 | 12.8 | clots < 3 min |
| legacy | 200 | 1.2 | 0.4 | 680.2 | 12.3 | clots < 3 min |
| legacy | 400 | 1.2 | 0.4 | 683.8 | 11.6 | clots < 3 min |

- model: threshold all clot (measured 50-100 um) -> FALLA
- legacy: threshold all clot (measured 50-100 um) -> FALLA
- j_X = 4.0e-13 mol m^-2 s^-1 from Gamma_TF = 8e-13 mol/m^2 and kappa_TF = 0.5 1/s (E1), not fitted.

References:
- Kastrup CJ, Runyon MK, Shen F, Ismagilov RF. Modular chemical mechanism predicts spatiotemporal dynamics of initiation in the complex network of hemostasis. PNAS 103 (2006) 15747-15752
- Kastrup CJ, Shen F, Runyon MK, Ismagilov RF. Characterization of the threshold response of initiation of blood clotting to stimulus patch size. Biophys J 93 (2007) 2969-2977
