# LeftAtrium2DPTAGX

`LeftAtrium2DPTAG` con tres interruptores nuevos (todos activos por defecto):

| Interruptor | Opción | Qué cambia |
|---|---|---|
| Factor Xa | `-ptag_xa true` | La iniciación pasa por $X$ transportado: $\mathcal L X=-k_{dX}X$, emisión $j_X E$ desde el endocardio (en la capa $\delta$ como fuente $(j_X/\delta)\chi\tilde E$, o como flujo de Neumann si $\delta=0$), y $R=(k_aT^2/(K_A+T)+\kappa X)P/P_{ref}$. El canal legado $P\to T$ de pared se apaga. |
| Transporte antisimétrico | `-ptag_skew true` | Añade $\tfrac12\int(\nabla\cdot\mathbf u_h)(c-c_{in})v$ al transporte de cada especie. Corrige una inestabilidad del driver original (ver abajo). |
| Fibrina inmóvil | `-ptag_immobile_fibrin true` | $F$ es un campo nodal, $\partial_tF=R_F$, en vez de $G_0-G$ (que el flujo advecta). El Brinkman lee este $F$. $G+F$ se conserva en el paso nodal. |

Con `-ptag_xa false -ptag_immobile_fibrin false -ptag_skew false` el driver es `LeftAtrium2DPTAG`: verificado bit a bit (las 17 columnas del CSV idénticas en 240 pasos, SR, malla gruesa).

## Constantes nuevas

| Opción | Símbolo | Defecto | Nota |
|---|---|---|---|
| `-ptag_j_x` | $j_X$ (mol m⁻² s⁻¹) | $j_0k_{dX}/\kappa=6.7\times10^{-15}$ | Valor "igualado": en estasis $X\to j_X\chi\tilde E/(\delta k_{dX})$ y la iniciación coincide con la legada; solo difieren donde el flujo lava $X$. Si se cambia `-ptag_wall_flux`, el igualado lo sigue. |
| `-ptag_kappa_x` | $\kappa$ (1/s) | 3 | Conversión por protrombinasa. Entra solo como $\kappa j_X$. |
| `-ptag_k_dx` | $k_{dX}$ (1/s) | $10^{-2}$ | Inhibición plasmática de Xa; longitud de lavado $\sqrt{D_P/k_{dX}}\approx68$ µm, tiempo $1/k_{dX}=100$ s. |

También quedan expuestas: `-ptag_k_on`, `-ptag_P0`, `-ptag_A0`, `-ptag_tau_a` (umbral de activación endotelial, Pa), `-ptag_e_osi` (`false`: $E=\sigma(\mathrm{TAWSS})$ sin el factor $\min(2\,\mathrm{OSI},1)$), `-ptag_csv`, `-ptag_xdmf` y `-ptag_laa_box xmin,xmax,ymin,ymax`.

## Salidas

CSV: las 17 columnas de `LeftAtrium2DPTAG` y además
`maxX, intX, intT, intF, laaX, laaT, laaF` (integrales por unidad de profundidad, mol/m; las `laa*` sobre la caja `-ptag_laa_box`, cero si no se da).
XDMF: los mismos campos y `factorXa`. `fibrin` es el $F$ que usa el Brinkman (inmóvil o $G_0-G$ según el interruptor).

La caja del apéndice se saca de la malla (bounding box de las aristas 2, 3, 4):

```bash
python3 tools/laa_box.py resources/LA2D_coarse.mesh          # -> -ptag_laa_box -0.0517,-0.0277,-0.0039,0.0039
```

## Compilar

`CMakeLists.txt` de `examples/Heart/Thrombosis` ya incluye `add_subdirectory(LeftAtrium2DPTAGX)` (bajo PETSc+MPI).

```bash
cd build && cmake . && make LeftAtrium2DPTAGX
```

## Casos para comparar

Cada caso en su directorio (o con `-ptag_csv`/`-ptag_xdmf` distintos). Desde `build/`; `X` es el ejecutable, `S` esta carpeta (rutas absolutas, porque cada caso corre en su directorio).

```bash
X=$PWD/examples/Heart/Thrombosis/LeftAtrium2DPTAGX/LeftAtrium2DPTAGX
S=$(cd ../examples/Heart/Thrombosis/LeftAtrium2DPTAGX && pwd); R=$S/resources
M=$R/LA2D_coarse.mesh                      # o LA2D_idealised.mesh (h = 1.2 mm)
BOX=$(python3 $S/tools/laa_box.py $M)
SR="-ptag_pv $R/presion_PV_SR.dat -ptag_mv $R/presion_MV_SR.dat"
AF="-ptag_pv $R/presion_PV_AF.dat -ptag_mv $R/presion_MV_AF.dat"
O="-ptag_mesh $M $BOX -ptag_dt 1e-2 -ptag_flow_cycles 2 -ptag_species_cycles 150 -ptag_replay true -ptag_species_substeps 2 -ptag_output_every 0"

c() { mkdir -p $1 && (cd $1 && shift && mpirun -n 4 $X $O "$@" > log.txt 2>&1); }

# 1. Ritmo: SR vs FA con el modelo nuevo
c sr_x $SR ;  c af_x $AF
# 2. Qué aporta cada interruptor (en FA)
c af_leg   $AF -ptag_xa false -ptag_immobile_fibrin false   # cinética de LeftAtrium2DPTAG, transporte estable
c af_orig  $AF -ptag_xa false -ptag_immobile_fibrin false -ptag_skew false   # = LeftAtrium2DPTAG
c af_xonly $AF -ptag_immobile_fibrin false                    # Xa, F advectada
c af_fonly $AF -ptag_xa false                                 # sin Xa, F inmóvil
# 3. Dosis de factor tisular / endotelio
c af_jx3   $AF -ptag_j_x 2e-14
c af_tau   $AF -ptag_tau_a 0.2        # endotelio más tolerante al cizalle bajo
c af_eosi  $AF -ptag_e_osi false      # E sin OSI: la estasis unidireccional también activa
# 4. Anticoagulantes (escalas ilustrativas, no farmacocinética)
c af_apix  $AF -ptag_kappa_x 0.6                    # anti-Xa directo (apixabán/rivaroxabán): kappa x0.2
c af_dabi  $AF -ptag_k_a 0.02                       # anti-IIa directo (dabigatrán): amplificación x0.2
c af_hep   $AF -ptag_k_on 35.5 -ptag_k_dx 5e-2      # heparina: AT x5 sobre T y sobre Xa
c af_warf  $AF -ptag_P0 0.42e-3 -ptag_j_x 2e-15     # warfarina: II y X al 30 % (P_ref fijo)

python3 $S/tools/compare_runs.py sr_x af_x af_leg af_xonly af_fonly af_apix af_dabi af_hep af_warf --png compare.png
```

`compare_runs.py` imprime, por caso: $\max T$, $\int\!\!\int T$ (exposición a trombina), la misma integral en el apéndice, $\max X$, la fibrina total y la del apéndice al final, su cociente y la deriva máxima de $P+T-A$. Con `--png` grafica las integrales en el tiempo. Acepta también CSV de `LeftAtrium2DPTAG` (deja vacías las columnas que no tiene).

## Verificación (contenedor, malla gruesa, 1 proceso)

- Modo legado (`-ptag_xa false -ptag_immobile_fibrin false -ptag_skew false`) = `LeftAtrium2DPTAG`, bit a bit (240 pasos, 17 columnas).
- Paso nodal aislado: $P+T-A$ y $G+F$ se conservan a $10^{-16}$, todo $\ge0$ incluso con $\Delta t=50$ s. En estasis, el $j_X$ igualado llega a la misma meseta de $T$ que el canal legado, con un retraso de $1/k_{dX}$.
- Corridas acopladas: $\Delta t=10^{-2}$, 2 + 40 ciclos (33.6 s), replay, `substeps 2`, cinética acelerada (`-ptag_k_a 0.3 -ptag_wall_flux 2e-10`), caja LAA. Todas con `-ptag_skew true`; la deriva queda en $\le1.8\times10^{-4}$.

| caso | max T (nM) | ∫∫T | LAA ∫∫T | max X (pM) | ∫F final | LAA ∫F final |
|---|---|---|---|---|---|---|
| SR, Xa+F inmóvil | 1.49 | 4.0e-10 | 2.9e-10 | 28 | 1.1e-8 | 8.1e-9 |
| FA, Xa+F inmóvil | 1.51 | 2.4e-10 | 1.3e-10 | 29 | 6.9e-9 | 3.5e-9 |
| SR, idem, $k_{dX}=0.1$ | 34 | 3.8e-9 | 2.6e-9 | 93 | 9.9e-8 | 6.9e-8 |
| FA, idem, $k_{dX}=0.1$ | 35 | 2.0e-9 | 1.0e-9 | 94 | 5.3e-8 | 2.7e-8 |
| SR, idem, $k_{dX}=0.1$, E sin OSI | 39 | 8.8e-9 | 6.3e-9 | 136 | 2.2e-7 | 1.5e-7 |
| FA, idem, $k_{dX}=0.1$, E sin OSI | 38 | 1.0e-8 | 6.6e-9 | 191 | 2.5e-7 | 1.6e-7 |
| SR, cinética legada | 203 | 1.4e-8 | 9.5e-9 | – | 2.2e-7 | 1.7e-7 |
| FA, cinética legada | 200 | 6.4e-9 | 3.5e-9 | – | 1.2e-7 | 6.8e-8 |

Lectura:
1. Con $k_{dX}=10^{-2}$, en 33 s $X$ llega solo al ~27 % de su meseta; el canal legado es instantáneo. Por eso Xa da 1.5 nM frente a 200 nM. Para comparar regímenes hay que correr $\gtrsim150$ ciclos.
2. Con $k_{dX}=0.1$ ($X$ en régimen en ~10 s), Xa produce ~45 % de la fibrina legada en ambos ritmos: es el lavado de $X$.
3. **El orden SR/FA lo decide $E$, no Xa.** Con $E=\sigma(\mathrm{TAWSS})\min(2\,\mathrm{OSI},1)$ el apéndice en SR tiene OSI medio 0.31 y en FA 0.14, así que SR activa 2.2 veces más pared del apéndice y forma más fibrina ahí (al revés de lo clínico). Con `-ptag_e_osi false` FA ≥ SR en exposición a trombina y $\max X$ (191 vs 136 pM: Xa se acumula más en la estasis), y el apéndice iguala. La forma de $E$ es la decisión de modelado más influyente en estos casos.

## Advertencias

- **Inestabilidad del transporte en `LeftAtrium2DPTAG`.** P1/P1 con VMS deja $\nabla\cdot\mathbf u_h\sim\pm10$ 1/s en los tramos cortos de pared entre dos venas (2–3 mm, uno o dos elementos). Con la forma $\int(\mathbf u\cdot\nabla c)v$ una desviación del estado de entrada crece ahí a tasa $\sim|\nabla\cdot\mathbf u_h|/2$. En SR (malla gruesa, $\Delta t=10^{-2}$) $P$ y $A$ bajan al 15 % en 10 ciclos en un nodo, $G$ al 12 %, y $G_0-G$ marca fibrina ($\max F=G_0$ en el ciclo 12). Con fibrina inmóvil y Xa sigue creciendo: $G$ llega a $10^3G_0$ en el ciclo 37. La forma antisimétrica pura no basta: estabiliza, pero fuerza el estado de entrada con $\tfrac12(\nabla\cdot\mathbf u_h)c_{in}$ ($P$ hasta $3.6P_0$). Hecha en torno a $c_{in}$, las constantes siguen siendo exactas y la desviación no se amplifica: la deriva queda acotada en $\approx1.2\times10^{-4}$ (9 % de $P_0$, sobreoscilaciones de SUPG). **Los resultados de especies en SR del driver original deben revisarse.**

- Con F advectada, `maxF` arranca en ~0.3 µM en cuanto empieza el transporte, aun sin trombina: son las oscilaciones de SUPG en $G$, que $G_0-G$ lee como fibrina (y que el Brinkman ve, con $K_F=2$ µM). La fibrina inmóvil solo crece con $R_F\ge0$.
- `invariantDrift` mide la deriva puntual de $P+T-A$; viene del transporte (SUPG no es positivo y el paso nodal recorta los negativos), no de la cinética, que conserva exactamente (verificado en el paso nodal aislado, también con $\Delta t=50$ s).
- Los parámetros de Xa y las escalas de los fármacos son de literatura u orden de magnitud; $j_X$ y $k_a$ siguen siendo las constantes por calibrar.
- La cinética de Xa tiene tiempo $1/k_{dX}=100$ s: hacen falta $\gtrsim150$ ciclos de especies (120 s) para que $X$ llegue a régimen en estasis.

## Opciones para benchmarks (Qureshi 2021, `../TestPTAGX/Q1_Qureshi2021`)

| Opción | Qué hace |
|---|---|
| `-ptag_inlet_velocity V` | Velocidad prescrita en las venas, u = −V(t) n, con la onda de `-ptag_pv` reescalada a pico V (m/s). Sin término de presión en las venas. |
| `-ptag_laa_wall_velocity W`, `-ptag_laa_wall_phase φ` | Pared de la orejuela (etiquetas 2, 3, 4) con velocidad normal −W sin(2πt/T + φ) (transpiración sobre malla fija). |
| `-ptag_wall_penalty β` | Penalización de las dos condiciones anteriores (por defecto 1e5 Pa·s/m). |
| `-ptag_seed x,y,r,T` | Trombina T (mol/m³) en un disco de radio r al inicio. |
| `-ptag_mu μ` | Sangre newtoniana. |
| `-ptag_k_F`, `-ptag_K_G`, `-ptag_D_P`, `-ptag_D_G` | Constantes de la conversión de fibrinógeno y difusividades. |
| `-ptag_laa_poly archivo` | Polígono de la orejuela (una línea "x y" por vértice) para las columnas `laa*`; reemplaza la caja. |
| `-ptag_csv_every k` | Fila del CSV cada k pasos (y al final de cada ciclo). |

Columnas nuevas del CSV: `laaArea, laaMeanU, laaMaxF, laaMaxT`. Con las opciones nuevas en sus valores por defecto, el modo legado sigue idéntico bit a bit a `LeftAtrium2DPTAG` (verificado de nuevo).

## Amplificación superficial (`-ptag_sigma_s`)

En la capa endocárdica (δ = `-ptag_layer`), θ = 1 + (σ/δ)·χ·Ẽ multiplica a k_a, con χ el mismo factor de dosis de la capa, así que ∫(θ − 1) dV = σ ∫E dA. σ (m) es la membrana procoagulante por unidad de área de endocardio activado. Por defecto σ = 0, idéntico al legado bit a bit. Validación y calibración: `../TestPTAGX` (σ = 2×10⁻³ m sale de Shen 2008; ver `RESULTS.md`, que muestra que ese valor rompe el umbral sin flujo de Kastrup).
