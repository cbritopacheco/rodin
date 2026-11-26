# Tests de validación del modelo de aurícula PTAGX

Cada carpeta es un test con su script, sus datos de referencia y su criterio. Un test escribe sus figuras y un `summary.md` en `results/`. Los resultados y su lectura están en `RESULTS.md`.

| Carpeta | Contra qué | Tipo de dato | Qué prueba del modelo | Costo |
|---|---|---|---|---|
| `E1_CAT` | Trombograma (CAT), plasma pobre en plaquetas, 5 pM FT, 4 µM fosfolípido; intervalos ELSA-Brasil 2022 | experimental | cinética 0D: lag, pico, ETP; qué k_a admite | segundos |
| `E2_Thrombodynamics` | Thrombodynamics (crecimiento del coágulo desde una superficie con FT, sin flujo): Tlag, Vi, Vst, CS30 | experimental | propagación espacial de la trombina (onda) | 7 × ~3 min |
| `E3_DOAC` | Thrombodynamics-4D en pacientes con FA bajo rivaroxabán, apixabán y dabigatrán, pico vs valle (ADAGE) | experimental (clínico) | el módulo Xa: anti-Xa frente a anti-IIa | 13 × ~3 min |
| `E4_PatchThreshold` | Umbral de tamaño de parche de FT para coagular plasma en reposo (Kastrup 2006): 50–100 µm | experimental | iniciación desde una lesión pequeña, sin ajustar nada (axisimétrico: parche = disco) | 10 × ~8 min |
| `E5_ShearThreshold` | Plasma con flujo sobre un parche de FT de 200 µm (Shen 2008): inicia a 5 s⁻¹, no a 40 s⁻¹ | experimental | iniciación contra el lavado | 16 × ~1 min |
| `Q1_Qureshi2021` | Aurícula 2D con cuatro orejuelas (CW, WS, BR, CA), SR y FA (Qureshi et al., CinC 2021, Tabla 1) | numérico | flujo y lavado en la orejuela; formación de trombo | calibración ~30 min; cada caso de 60 s, ~45 min |

Los tests E2–E5 usan `PTAGXBench`, una copia de `PTAGP` (en `../PTAGXBench`, sin tocar tu `PTAGP`) con la capa de amplificación superficial y una opción axisimétrica (`-pg_axisym`, para E4). Con las opciones de `ptagx_common.PTAGX` (plaquetas apagadas, θ₀ = 1, Xa activo, fibrina inmóvil sin monómero ni tope, sin trombomodulina), su paso nodal es el de `LeftAtrium2DPTAGX`. Q1 usa `LeftAtrium2DPTAGX` directamente.

**Qué se calibró y con qué.** k_a sale de E2. Todo lo demás es predicción: E1, E3, E4, E5 y Q1. (σ del juego `SURFACE` sale de E5 a 5 s⁻¹.)

**Veredictos.** Cada test da su veredicto sobre `MODEL`. Los otros juegos aparecen en el mismo `verdict.md` como sensibilidad, sin entrar en el veredicto.

## Cómo se corre

Compilar los dos ejecutables (desde la raíz del repo):

```bash
cd build && cmake . && make PTAGXBench LeftAtrium2DPTAGX
cd ../examples/Heart/Thrombosis/TestPTAGX
```

Todo, con el veredicto de cada test apenas termina y un resumen al final:

```bash
python3 run_all.py --jobs 8                       # E1-E5 y Q1 completos
python3 run_all.py --skip-q1 --jobs 8             # solo los experimentales (~1 h en 2 núcleos)
python3 run_all.py --only Q1 --jobs 8 --q1-kinetics qureshi,ptagx,ptagx_surface,ptagx_legacy
```

Un test suelto (cada uno imprime su veredicto al final y lo guarda en `results/verdict.md`):

```bash
python3 E1_CAT/run_E1.py
python3 E2_Thrombodynamics/run_E2.py --jobs 8
python3 E3_DOAC/run_E3.py --jobs 8
python3 E4_PatchThreshold/run_E4.py --jobs 8
python3 E5_ShearThreshold/run_E5.py --jobs 8
python3 E5_ShearThreshold/calibrate_jx.py         # dosis de FT que reproduce el umbral de cizalle
python3 E5_ShearThreshold/calibrate_sigma.py      # (histórico) recalibrar sigma de SURFACE
python3 E4_PatchThreshold/scan_sigma.py           # sigma que tolera el parche de 50 um sin flujo (tira 2D)
python3 Q1_Qureshi2021/run_Q1.py --stage calib --jobs 2
python3 Q1_Qureshi2021/run_Q1.py --stage coag --jobs 8 --kinetics qureshi,ptagx,ptagx_surface,ptagx_legacy
python3 Q1_Qureshi2021/run_Q1.py --stage report   # veredicto por cinética en results/<cinética>/verdict.md
```

Opciones comunes: `--jobs`, `--force` (repetir casos existentes), `--plot-only`, `--quick`. Variables de entorno: `PTAGP_EXE` (ejecutable `PTAGXBench`), `PTAGX_EXE` (`LeftAtrium2DPTAGX`), `PTAGP_MPIRUN` (en un contenedor: `"mpirun --allow-run-as-root -n 1"`).

**Parámetros que se prueban** (`ptagx_common.py`):
- `MODEL` (propuesto, el de los veredictos): k_a de E2, κ = 10 s⁻¹, ε = 0.1, σ = 0.
- `SURFACE` (sensibilidad histórica): κ = 3, ε = 1 y σ = 2.0×10⁻³ m (de E5, con `calibrate_sigma.py`). Arregla E5 y rompe E2–E4; reemplazado por la trombomodulina (ver `RESULTS.md`).
- `LEGACY`: los valores por defecto de PTAGX (k_a = 0.10, κ = 3, ε = 1, σ = 0).

## El modelo

Especies transportadas por L c = ∂_t c + u·∇c − ∇·(D ∇c): P (protrombina), T (trombina), A (antitrombina), G (fibrinógeno), X (**factor Xa**, no el zimógeno), f (monómero de fibrina). F (gel) es nodal.

| Reacción | Expresión |
|---|---|
| amplificación | R = [k_a θ T²/(K_A+T) + κ X θ_X] P/P_ref |
| cofactor Va | θ_X = ε + (1−ε) T/(K_A+T) |
| inhibición | R_I = k_on A T |
| corte del fibrinógeno | R_F = k_F T G/(K_G+G) → **monómero** |
| polimerización | R_pol = k_pol f (f+F) → gel |
| iniciación (pared) | −D ∂_n X = j_X E, con j_X = Γ_FT·k_cat |
| protección (pared) | k_TM (1−E) T, pagada por A |

- **X es Xa.** El flujo de pared j_X = Γ_FT·k_cat lleva el FT y el FVIIa; no hacen falta como campos. El FX del plasma (170 nM) pierde < 1 nM por corrida, así que llevarlo sería un campo y una K_M de más.
- **El cofactor θ_X** (ε, `-pg_eps_x` / `-ptag_eps_x`) es lo único nuevo en la cinética: Xa libre activa mal la protrombina, Xa·Va·fosfolípido ~10⁵ veces mejor, y el Va lo hace la trombina. Así la iniciación es un umbral y la propagación (θ_X → 1) no cambia. ε = 1 recupera el modelo anterior bit a bit.
- **El monómero** (`-pg_fibrin_monomer` / `-ptag_fibrin_monomer`, k_pol = 10³ m³ mol⁻¹ s⁻¹) permite que la trombina corte fibrinógeno lejos de la lesión sin dejar gel ahí. Sin flujo no cambia nada.
- **La trombomodulina** (`-ptag_k_tm`) es la otra cara de E: el mismo E que enciende el FT apaga la protección. La protección es la extensión de (1−E), o sea L̃ − Ẽ con L̃ la extensión de 1; (1−Ẽ) no es una capa, vale ~1 en toda la luz. `-ptag_tfpi_loss` hace la mitad del TFPI, k_dX → k_dX(1−λẼ); default 0.
- **Amplificación superficial (opción 2, en desuso).** θ = 1 + (σ/δ)·χ·w, ∫(θ−1) dV = σ∫E dA. En `LeftAtrium2DPTAGX` (`-ptag_sigma_s`) w = Ẽ; en `PTAGXBench` (`-pg_sigma_s`, `-pg_delta_s`, `-pg_patch`) w = exp(−d/δ_s). Con σ = 0 ambos son idénticos bit a bit a `LeftAtrium2DPTAG` y a `PTAGP`.

**Axisimétrico (`-pg_axisym`, solo sin flujo).** x es el radio r, eje en x = 0; toda integral de transporte lleva el peso r (dV = r dr dy). Un parche [0, a] en y = 0 es un disco de radio a. Verificado contra la solución exacta de un disco con flujo uniforme y decaimiento lineal en un semiespacio: 0.4 % (a = 25 µm) y 0.1 % (a = 100 µm), una vez descontado el error de separación de operadores (D·(1 + k dt)).

## Q1: qué es exacto y qué se supuso

- **Geometría.** WS y BR salen de los STL de los autores (github.com/ahmedqureshi1/2D_LAA); el cuerpo auricular es común a ambos a 0.15 mm. CW se digitalizó del panel principal de la Fig. 1 (registrado sobre el cuerpo, rms 0.42 mm). CA se digitalizó de su miniatura; la figura no da su escala, y se fijó para que el área de la orejuela sea la media de las otras tres (4.2 cm²). Los parches PV1–4, la mitral y el parche LAA, y las posiciones P1–P3, salen de la Fig. 1.
- **Cinética de Qureshi como límite de PTAGX.** Con D_G = D_Fn, la fibrina advectada F = G₀ − G es exactamente su ecuación de Fn. Con K_G = 1 mol/m³ ≫ G y k_F = K_eff K_G, la conversión es K_eff Fg Th. La trombina queda como trazador pasivo (k_a = k_on = 0, sin fuente de pared): su R_Th (modelo de Ataullakhanov) no está en el artículo. Este límite aísla el transporte (lavado frente a estasis).
- **Supuesto: forma de onda en las venas.** Ondas S, D y, en SR, reversión A; período 1 s; pico 7.5 cm/s (SR) y 5 cm/s (FA). En `waveforms/`.
- **Supuesto: pared de la orejuela.** Velocidad normal sinusoidal (transpiración sobre malla fija, sin ALE), con la contracción a t = 0.9 s. Su amplitud se calibra en CW para que la velocidad media en la orejuela sea la de la Tabla 1 (0.12 cm/s SR, 0.022 cm/s FA). Las otras tres formas predicen su V_mean.
- **Opciones nuevas de `LeftAtrium2DPTAGX`.**
  - Contorno: `-ptag_inlet_velocity` (velocidad en las venas) y `-ptag_laa_wall_velocity`, `-ptag_laa_wall_phase`, `-ptag_wall_penalty` (pared de la orejuela).
  - Condición inicial y cinética: `-ptag_seed x,y,r,T`, `-ptag_mu`, `-ptag_k_F`, `-ptag_K_G`, `-ptag_D_P`, `-ptag_D_G`.
  - Salida: `-ptag_laa_poly` (polígono de la orejuela) y `-ptag_csv_every`.
  - Columnas nuevas del CSV: `laaArea`, `laaMeanU`, `laaMaxF`, `laaMaxT`.

## Lo que no se incluyó, y por qué

- **Okorie 2008 y Colace 2012** (umbral de FT vs cizalle, fibrina y plaquetas sobre colágeno). Son en sangre entera y su respuesta depende de las plaquetas; PTAGX no las tiene. Están en `../Test` para PTAGP.
- **Taylor 2014 (escalón, MRI de trombo).** Los valores del tamaño del trombo no están en el resumen, y el trombo es mayormente plaquetario.
- **Velocidades de vaciado de la orejuela por ecocardiografía** (SR > 40 cm/s, FA < 20 cm/s). Una aurícula 2D idealizada no se compara en valor absoluto.
