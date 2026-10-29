# Válvulas de cuchillo FPC8 / FPC16 — fuerza del fluido sobre la compuerta

Navier–Stokes incompresible, P1/P1 estabilizado (OSS split, proyección
ortogonal de `examples/Heat`), Picard, viscosidad power-law. Mallas sólo de
fluido: tubería recta (diámetro interior de la válvula) menos el cuchillo
desplazado a la apertura pedida.

## Mallas (`resources/examples/Valve/`)

| archivo | D [mm] | apertura | gap/D | vértices | tetraedros | γ mín |
|---|---|---|---|---|---|---|
| FPC8_open000.mesh  | 203.2 | 0 % (asentada)   | −0.017 | 86 990 | 435 861 | 0.21 |
| FPC8_open025.mesh  | 203.2 | 25 %             | 0.236 | 85 581 | 431 394 | 0.14 |
| FPC8_open050.mesh  | 203.2 | 50 %             | 0.488 | 77 431 | 387 755 | 0.12 |
| FPC8_open075.mesh  | 203.2 | 75 %             | 0.740 | 73 793 | 365 277 | 0.13 |
| FPC8_open100.mesh  | 203.2 | 100 % (retirada) | 0.992 | 77 560 | 374 547 | 0.11 |
| FPC16_open000.mesh | 382.0 | 0 %              | −0.013 | 93 696 | 454 547 | 0.20 |
| FPC16_open025.mesh | 382.0 | 25 %             | 0.238 | 87 128 | 429 216 | 0.15 |
| FPC16_open050.mesh | 382.0 | 50 %             | 0.490 | 78 514 | 386 738 | 0.11 |
| FPC16_open075.mesh | 382.0 | 75 %             | 0.741 | 73 984 | 361 570 | 0.08 |
| FPC16_open100.mesh | 382.0 | 100 %            | 0.992 | 77 326 | 371 110 | 0.10 |

- **Dominio de fluido:** tubería recta del diámetro interior de la válvula
  **unida a la cavidad interior del revestimiento** (FPC8, `REVESTIMIENTO`)
  o de la camisa (FPC16, `CAMISA PU`) tomada del STEP, menos el cuchillo
  desplazado a la apertura. En el plano del cuchillo la cavidad tiene:
  - **abajo:** la ranura del asiento (radio R_k > R: 3.4 mm más profunda
    que el diámetro en el FPC8, 5 mm en el FPC16);
  - **arriba:** la ranura por la que sube el cuchillo, hasta 1 mm bajo el
    tope del revestimiento (ahí sella la prensaestopa). El cuchillo la
    llena por encima de su borde; quedan mojados los dos bolsillos
    laterales entre el diámetro y las paredes de la ranura bajo el borde
    del cuchillo y, en el FPC16, las cuñas entre el bisel y la pared de la
    ranura en toda su altura.
  - Cada bolsillo lateral termina en una cúspide (el costado del cuchillo
    es tangente a su borde y coincide con la pared): se recorta donde es
    más angosto que 0.0075 D (1.5 mm FPC8, 2.9 mm FPC16).
  - `--cavity groove` deja sólo la ranura del asiento (como si los labios
    de la camisa sellaran sobre el cuchillo); `--cavity none` vuelve a la
    tubería recta.
- **Consecuencia en la fuerza:** la presión de línea (150 psi) sobre las
  cuñas del bisel del FPC16 no se compensa con el otro lado del cuchillo:
  al 25 % da ≈ −7 kN en z, y la del borde ≈ +6 kN en y, aun sin flujo
  (malla gruesa de prueba). Es la parte `p_out ∫ n dS` que `F_dyn` resta;
  la parte del flujo es `F_dyn`.
- Unidades: **metros**. Flujo en +z, carrera del cuchillo en +y (+y abre).
- Tubería: 3D aguas arriba, 8D aguas abajo del plano del cuchillo.
- Apertura `f`: centro del borde semicircular del cuchillo en
  `y_c = f (R + R_k − δ)`, con `δ = 0.0079 D` la protrusión residual que
  trae el CAD del FPC8 abierto (1.6 mm). `f = 0` es el cuchillo asentado en
  la ranura del asiento: la malla tiene dos regiones desconectadas y el
  caso es puramente hidrostático (`F_z = Δp · A_mojada`).
- Tamaños: `h = D/80` en el cuchillo, `D/133` en las paredes de la cavidad
  y `D/160` en la unión cuchillo–pared, `D/26` en la caja del chorro
  (−0.75D … 3D aguas abajo), `D/10` lejos.
- Caras del cuchillo: se clasifican sobre una malla superficial gruesa
  (centro de cada triángulo desplazado según su normal, dentro o fuera del
  cuchillo del CAD), no por coordenadas.
- Atributos de frontera: 1 entrada, 2 salida, 3 pared (tubería y cavidad),
  4 cuchillo. Volumen: 1.
- Regenerar / cambiar resolución: `python3 make_meshes.py FPC8.stp FPC16.stp
  --hk 0.0125 --hm 0.0385 --hf 0.1 [--cavity full|groove|none]` (requiere
  `pip install gmsh`). Malla gruesa de prueba: `--hk 0.03 --hm 0.08 --hf 0.2
  --ldown 4` (~55 k tets).

## Ejecutable `KnifeGateValve`

```sh
cmake --build build -j --target KnifeGateValve
cd scratch
# compuerta abierta (25/50/75/100 %): agua a 3 m/s en una línea de 150 psi
mpirun -n 8 ../build/examples/Valve/KnifeGateValve \
    -kv_mesh ../resources/examples/Valve/FPC8_open050.mesh -kv_out fpc8_50 -kv_umean 3
# compuerta cerrada (0 %): 150 psi aguas arriba, 0 aguas abajo (hidrostático)
mpirun -n 8 ../build/examples/Valve/KnifeGateValve \
    -kv_mesh ../resources/examples/Valve/FPC8_open000.mesh -kv_out fpc8_0 \
    -kv_inlet pressure -kv_pin 1.0342e6 -kv_pout 0 -kv_tend 0.1
```

Valores por defecto: agua a 20 °C (ρ = 998.2 kg/m³, μ = 1.002e-3 Pa·s,
es decir `K = μ`, `n = 1`), U = 3 m/s, p_out = 150 psi = 1.0342 MPa.

Opciones (prefijo `-kv_`): `mesh`, `out`, `model` (`powerlaw`|`bingham`),
`rho`, `K` (Pa·sⁿ), `n`, `gdmin`, `mumin`, `mumax`, `tau_y` (Pa), `mu_p`
(Pa·s), `papa_m` (s), `bingham_ratio`, `inlet` (`velocity`|`pressure`),
`umean`, `q` (caudal en m³/s; reemplaza a `umean`), `profile`
(`auto`|`laminar`|`turbulent`|`uniform`), `pin`, `pout`, `ramp`,
`total_pressure`, `dt`, `dtmin`, `dtmax`, `tend`, `adaptive_dt`, `easy`,
`hard`, `maxit`, `tol`, `tau_dt`, `bdf`, `extrapolate`, `lumped`, `refactor_every`, `refactor_its`, `picard_monitor`, `tavg`, `steady_tol`,
`steady_window`, `output_every`.
Solver: `-kv_flow_*` (MUMPS LU por defecto; si no hay MUMPS,
`-kv_flow_pc_factor_mat_solver_type superlu_dist` o `-kv_flow_ksp_type gmres
-kv_flow_pc_type ilu` en secuencial), proyecciones `-kv_proj_*` (CG+Jacobi).

### Condiciones de borde

- **Entrada de velocidad (por defecto, compuerta abierta).** Perfil
  desarrollado con caudal exacto `U·A_in` (se normaliza con el interpolante
  P1 que impone el Dirichlet). `auto` elige con el Reynolds de
  Metzner–Reed: laminar si Re_MR < 2100 (perfil power-law
  `1 − (r/R)^((n+1)/n)`, Poiseuille si n = 1), si no ley 1/7. Agua a 3 m/s
  da Re ≈ 6e5 (FPC8) y 1.1e6 (FPC16): perfil turbulento. Rampa C¹ de la
  velocidad en `ramp` segundos.
- **Entrada de presión** (`-kv_inlet pressure`): `σn = −p_in n` con
  `u × n = 0` por penalización. Es la única opción para la compuerta
  **cerrada**, porque ahí no se puede imponer caudal.
- **Salida libre**: `σn = −p_out n`, *directional do-nothing* contra el
  reflujo. No-slip en pared y cuchillo.

### Qué significan los 150 psi

CWP 150 psi es la presión de trabajo de la línea (lo que aguanta el
cuerpo), no la caída de presión a través de la válvula abierta. Por eso:

- Abierta: se impone el caudal y la línea queda a 150 psi en la salida; la
  caída Δp y la fuerza salen del cálculo. Imponer 150 psi de diferencia
  con la válvula abierta daría velocidades irreales (decenas de m/s).
- Cerrada: los 150 psi sí son la diferencia que sostiene la compuerta:
  `F_z = Δp · A_mojada` ≈ 33.5 kN (FPC8) y 118.5 kN (FPC16).
- Una presión uniforme p_out sobre el cuchillo mojado no da fuerza en z
  (las dos caras se compensan) pero sí en y: el borde está mojado y los
  costados (dentro de la ranura de la camisa) no. Esa parte depende de la
  presión en la cámara de la válvula, que no está modelada. Por eso el
  CSV trae también **F_dyn = F − p_out ∫ n dS**, la fuerza hidrodinámica.
  En FPC8 al 50 % la parte de presión uniforme es ~1.85 kN en y, frente a
  unas decenas de N de F_dyn.

### Salida

`<out>.csv` (t, dt, iteraciones Picard, Q_in, Q_out, p medias, **F**
consistente, **F_direct**, **F_pressure**, corrección de unión, **F_dyn**,
ζ) y `<out>.xdmf` (velocidad, presión, viscosidad y la carga nodal sobre
las paredes `wall_load`, vector en Pa). Al final imprime el promedio
temporal de F y F_dyn sobre la segunda mitad del intervalo (`-kv_tavg`).
Con `-kv_steady_tol/-kv_steady_window` se detiene si (F, Q) se congelan.

### Paso de tiempo y costo

A Re ~ 1e6 el chorro bajo el cuchillo es inestable: la fuerza oscila y lo
que vale es el promedio. Picard deja de contraer cuando el CFL en el chorro
pasa de ~5 (en la malla gruesa de prueba, dt = 5e-3 divergía y 1.25e-3
convergía en 5 iteraciones); el programa reintenta con dt/2 y adapta dt
según las iteraciones. Valores por defecto: dt = 5e-4 s, tolerancia de
Picard 1e-3, t_final = 0.6 s, promedio desde 0.3 s. Son ~1000–2000 pasos;
con 200k–350k incógnitas y MUMPS, del orden de un día por caso en un
equipo de escritorio con 8 núcleos (estimación, no medido). Para una
primera pasada rápida, generar mallas más gruesas con
`make_meshes.py --hk 0.025 --hm 0.06 --hf 0.15`.

## Pulpa como plástico de Bingham (FPC8, `run_valves_bingham.sh`)

Valores máximos de la tabla de la pulpa (fila de 8"): ρ = 1770 kg/m³,
τ_y = 7.0 Pa, μ_p = 0.0182 Pa·s, Q = 188.8 m³/h (U = Q/A = 1.62 m/s).
Salida a presión atmosférica (p_out = 0 manométrica).

```sh
cmake --build build -j --target KnifeGateValve
mkdir -p ~/valve_runs && cd ~/valve_runs
bash ~/rodin/examples/Valve/run_valves_bingham.sh            # 32 procesos, un caso a la vez
NJOBS=3 bash ~/rodin/examples/Valve/run_valves_bingham.sh    # o tres casos a la vez, ~10 procesos c/u
```

Deja `runs_bingham/FPC8_open{000,025,075,100}/` (log, CSV, XDMF cada 10
pasos) y `runs_bingham/forces_summary.csv`. Se puede relanzar: los casos
terminados se saltan. No para en estado estacionario (`-kv_steady_tol 0`):
corre hasta 20 D/U = 2.5 s y promedia desde 8 D/U = 1.0 s. El
post-proceso es el mismo: `python3 postprocess_valves.py runs_bingham`.

- **Ley constitutiva.** Bingham regularizado (Papanastasiou):
  μ(γ̇) = μ_p + τ_y (1 − e^(−m γ̇))/γ̇, que tiende a μ_p + τ_y m en las
  zonas que no fluyen. Por defecto m = (1000 − 1) μ_p/τ_y = 2.6 s, es decir
  μ(0) = 1000 μ_p = 18.2 Pa·s: la tensión de fluencia se recupera a más de
  99 % para γ̇ > 2 s⁻¹, y la cizalla en el flujo es de cientos de s⁻¹.
  Más rígido: `-kv_papa_m 26` (μ(0) = 182 Pa·s; Picard más lento) o
  `-kv_bingham_ratio`.
- **Régimen.** Re_B = ρUD/μ_p ≈ 3.2e4, He = ρτ_yD²/μ_p² ≈ 1.6e6,
  Bn = τ_yD/(μ_pU) ≈ 48. El Reynolds crítico de Hanks para ese He es
  ≈ 1.8e4: el flujo en la tubería es turbulento (débilmente) y el perfil
  de entrada `auto` es la ley 1/7. En laminar (`-kv_profile laminar`)
  sería el de Buckingham–Reiner, con tapón de radio 0.73 R. El log imprime
  estos números al inicio.
- **0 % (cerrada).** No pasa caudal; se impone 150 psi (CWP) aguas arriba y
  0 aguas abajo (`PIN_CLOSED`, `POUT_CLOSED`). El fluido está en reposo y
  la fluencia no interviene: F_z = Δp · A_mojada.
- **Mallas** `resources/examples/Valve/wake/` (las mismas de arriba con la
  caja del chorro a `D/31` en vez de `D/26`, y hasta 4.5 D aguas abajo en
  vez de 3 D: `make_meshes.py --only FPC8 --hm 0.032 --wake-end 4.5
  --outdir wake`):

  | archivo | vértices | tetraedros | γ mín |
  |---|---|---|---|
  | FPC8_open000.mesh | 137 248 | 721 860 | 0.20 |
  | FPC8_open025.mesh | 135 874 | 717 370 | 0.16 |
  | FPC8_open075.mesh | 125 841 | 661 257 | 0.18 |
  | FPC8_open100.mesh | 130 041 | 673 435 | 0.15 |

  ~1.7 veces las celdas de las mallas anteriores (~545 k incógnitas).
- **Solver iterativo (por defecto en este script).** El LU de MUMPS de estas
  mallas pide 25–30 GB (crece como N^(4/3); el MUMPS de Ubuntu sólo trae el
  ordenamiento PORD). `-kv_solver iterative` usa FGMRES con un
  precondicionador de complemento de Schur por bloques [u | p]
  (`fieldsplit`, factorización `lower`, S ≈ A11 − A10 diag(A00)⁻¹ A01) y un
  ciclo V de BoomerAMG (hypre) en cada bloque; el precondicionador se
  reconstruye cada 20 resoluciones o cuando FGMRES pasa de 150 iteraciones.
  Medido en una malla de 187 k celdas (1 proceso): 2.0 GB frente a 3.3 GB
  con LU, y las mismas fuerzas (a 6 cifras) en los primeros pasos. Para
  las mallas de la estela se estiman ~10–12 GB en total. `SOLVER=direct`
  vuelve al LU. El log imprime la memoria por paso (`memory [GB]`).

## Cómo se calcula la fuerza (y por qué así)

`F = −∫_cuchillo σn dS`, fuerza del fluido sobre el cuchillo, ejes globales.

1. **Reacción consistente (la que hay que usar).** El residuo discreto de
   momento sin eliminar Dirichlet, `R_k(φ_i) = a((u_h,p_h), φ_i e_k)` con
   todos los términos de volumen (estabilización incluida), es igual a
   `∫_∂Ω (σn)_k φ_i`. Entonces `F_k = −Σ_{i∈cuchillo} (R_k(φ_i) + ∫_pared p n_k φ_i)`;
   el segundo término quita, en los nodos de la unión cuchillo–pared, la
   presión de las caras de pared que comparten el nodo (exacto con p en P1).
   No se evalúa ningún gradiente de velocidad en una cara de frontera.
   Verificación: en el residuo, los nodos interiores quedan en 1e-12 (la
   forma lineal reproduce exactamente la bilineal), y en la válvula cerrada
   `F_z = Δp · A_mojada` a 1e-6 relativo.
2. **Integral directa de superficie (control).** `∫ (p n − μ(∇u+∇uᵀ)n) dS`
   con `∇u` tomado del tetraedro adyacente a cada cara, evaluando `Jacobian(u)`
   por la ruta `Point` (que resuelve la cara a su celda incidente).

**Lo que no funciona:** `Jacobian(u)` dentro de un `BoundaryIntegral` se
evalúa con el elemento de la *cara* (`IntegrationPoint`), es decir, sólo
el gradiente tangencial. Con no-slip ese gradiente es idénticamente cero,
así que la parte viscosa de la fuerza sale 0 y sólo queda la presión. La
diferencia entre (1) y (2) es O(h) en la parte viscosa (gradiente P1 de
primer orden en la pared); ambas tienen que converger al mismo valor al
refinar — eso es lo que se controló en `conv_*` (ver abajo).

## Estabilización

Split OSS (Codina): `τ_M ρ((∇u)a − Π[(∇u^n)a^n], (∇v)a)`,
`(τ_M/ρ)(∇p − Π[∇p^n], ∇q)`, `τ_C(div u − Π[div u^n], div v)`, con
`τ_M = [(2/dt)² + (2|u^k|/h)² + (4ν^k/h²)²]^{−1/2}`,
`τ_C = μ^k + ρ|u^k|h/2` y `a` la media por celda de `u^k` (P0; la
convección de Galerkin mantiene `(∇u)u^k`). El término 2/dt va por
defecto: sin él, en agua a baja velocidad `τ_M ~ h²/4ν` (≈ 1.5 s con
h = 2.5 mm) supera a dt en órdenes de magnitud, las proyecciones
retrasadas penalizan con peso τ/dt ≫ 1 todo cambio dentro del paso, la
presión se atrasa respecto del flujo y Picard se estanca
(`-kv_tau_dt false` lo quita). Las proyecciones Π van sobre P1 con masa
concentrada (`-kv_lumped false`: L2 consistente) y se actualizan una vez
por paso (en el iterado hacen que Picard contraiga sólo como τ/(τ+dt)).
Coeficientes (μ, τ, a) en P0, exactos para u en P1.

Tiempo: BDF2 de paso variable (BDF1 en el primer paso; `-kv_bdf 1` para
Euler implícito). Picard parte de la extrapolación lineal de
(u^{n−1}, u^n) (`-kv_extrapolate false` la quita).

## Rendimiento

El log imprime en cada paso el tiempo por fase (`time [s]: setup,
coefficients, assemble, solve, forces, projections, output`) y al inicio el
balance de la partición (`[partition] ... max/mean`).

- **Partición.** Sin Scotch, el particionador de respaldo puede dejar casi
  toda la malla en el rank 0 (medido aquí: 41 510 de 41 514 celdas con 2 y
  3 procesos), y los demás procesos quedan ociosos. Con
  `-DRODIN_USE_SCOTCH=ON` queda balanceada (max/mean = 1.01). Si el log
  muestra max/mean > 1.5, ese es el problema.
- **Ensamblaje.** Es ~80 % del paso: cada integrador recorre la malla por
  separado (~130 µs por celda y ensamblaje en 2 procesos). La factorización
  MUMPS es ~10 % en la malla gruesa de prueba.
- **Solver lineal.** GMRES precondicionado con factores MUMPS que se
  reutilizan entre iteraciones de Picard y entre pasos; se refactoriza cada
  20 resoluciones (`-kv_refactor_every`) o cuando GMRES pasa de 15
  iteraciones (`-kv_refactor_its`). En FPC16 25 % grueso: 1.9 → 0.4 s de
  resolución por paso, fuerzas idénticas a 1e-6. `-kv_refactor_every 1
  -kv_flow_ksp_type preonly` es la resolución directa de antes. El log dice
  cuántas factorizaciones e iteraciones de GMRES hubo en cada paso.
- Medido en la malla gruesa de prueba (FPC8 25 %, 41 k tetraedros, 2
  procesos, Picard = 2): 16.2 s por paso antes, 11.5 s con advección P0 en
  la OSS, reacción en un solo ensamblaje y proyección concentrada, 6.7 s
  además con la partición de Scotch.

## Notas de uso

- Lee el **promedio** temporal: a Re alto la estela es inestable y la
  fuerza oscila.
- Con la válvula cerrada usa `-kv_inlet pressure`; la salida es un recinto
  con `p = p_out` constante y ζ no se reporta (Q ≈ 0).
- Convergencia de malla (FPC8 50 %, dominio corto 1D + 2D, K = 10 Pa·s,
  ρ = 1000, Δp = 2 kPa, Re ~ 1: el caso más desfavorable, porque la parte
  viscosa pesa ~15 % de F_z). Estado estacionario:

  | h en el cuchillo | F_y reacción | F_y directa | F_z reacción | F_z directa | F_z sólo presión |
  |---|---|---|---|---|---|
  | D/25 | 0.589 | 1.168 | 21.36 | 19.45 | 17.71 |
  | D/40 | 0.473 | 0.885 | 20.82 | 19.55 | 17.81 |
  | extrapolado (orden 1) | 0.28 | 0.41 | 19.9 | 19.7 | — |

  Las dos evaluaciones convergen al mismo valor (la brecha baja como h:
  1.91 → 1.27 N en z) y el extrapolado coincide a 1 % en F_z. La fuerza
  lateral F_y es chica en este caso (~2 % de |F|). Las mallas de producción
  usan D/80 en el cuchillo. A Re alto la presión domina y la diferencia
  entre ambas se achica; si las dos discrepan más de unos pocos %, la malla
  del cuchillo no está resuelta para ese caso.
