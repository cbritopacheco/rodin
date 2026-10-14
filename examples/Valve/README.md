# Válvulas de cuchillo FPC8 / FPC16 — fuerza del fluido sobre la compuerta

Navier–Stokes incompresible, P1/P1 estabilizado (OSS split, proyección
ortogonal de `examples/Heat`), Picard, viscosidad power-law. Mallas sólo de
fluido: tubería recta (diámetro interior de la válvula) menos el cuchillo
desplazado a la apertura pedida.

## Mallas (`resources/examples/Valve/`)

| archivo | D [mm] | apertura | gap/D | vértices | tetraedros |
|---|---|---|---|---|---|
| FPC8_open000.mesh  | 203.2 | 0 % (asentada)   | −0.017 | 86 839 | 434 981 |
| FPC8_open025.mesh  | 203.2 | 25 %             | 0.236 | 78 224 | 402 177 |
| FPC8_open050.mesh  | 203.2 | 50 %             | 0.488 | 67 814 | 352 243 |
| FPC8_open075.mesh  | 203.2 | 75 %             | 0.740 | 58 601 | 306 189 |
| FPC8_open100.mesh  | 203.2 | 100 % (retirada) | 0.992 | 50 681 | 268 639 |
| FPC16_open000.mesh | 382.0 | 0 %              | −0.013 | 87 216 | 436 337 |
| FPC16_open025.mesh | 382.0 | 25 %             | 0.238 | 77 339 | 397 989 |
| FPC16_open050.mesh | 382.0 | 50 %             | 0.490 | 67 062 | 347 980 |
| FPC16_open075.mesh | 382.0 | 75 %             | 0.741 | 57 979 | 303 528 |
| FPC16_open100.mesh | 382.0 | 100 %            | 0.992 | 50 415 | 267 521 |

- Unidades: **metros**. Flujo en +z, carrera del cuchillo en +y (+y abre).
- Tubería: 3D aguas arriba, 8D aguas abajo del plano del cuchillo.
- Apertura `f`: centro del borde semicircular del cuchillo en
  `y_c = f (R + R_k − δ)`, con `δ = 0.0079 D` la protrusión residual que
  trae el CAD del FPC8 abierto (1.6 mm). `f = 0` es el cuchillo asentado en
  la ranura del asiento (el gap es negativo porque el borde entra en la
  ranura: la malla tiene dos regiones desconectadas y el caso es puramente
  hidrostático, útil para verificar `F_z = Δp · A_mojada`).
- Tamaños: `h = D/80` en el cuchillo (`D/160` en la unión cuchillo–pared,
  donde el bisel del FPC16 deja cuñas delgadas), `D/26` en la caja del chorro
  (−0.75D … 3D aguas abajo), `D/10` lejos. Calidad γ mínima ≥ 0.20 en todas.
- Atributos de frontera: 1 entrada, 2 salida, 3 pared de la tubería,
  4 cuchillo. Volumen: 1.
- Regenerar / cambiar resolución: `python3 make_meshes.py FPC8.stp FPC16.stp
  --hk 0.0125 --hm 0.0385 --hf 0.1` (requiere `pip install gmsh`). Para una
  malla más gruesa de prueba `--hk 0.04 --hm 0.08 --hf 0.2 --ldown 5`
  (~45 k tets, corre en segundos por paso en 1 núcleo).

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

Opciones (prefijo `-kv_`): `mesh`, `out`, `rho`, `K` (Pa·sⁿ), `n`, `gdmin`,
`mumin`, `mumax`, `inlet` (`velocity`|`pressure`), `umean`, `profile`
(`auto`|`laminar`|`turbulent`|`uniform`), `pin`, `pout`, `ramp`,
`total_pressure`, `dt`, `dtmin`, `dtmax`, `tend`, `adaptive_dt`, `easy`,
`hard`, `maxit`, `tol`, `tau_dt`, `picard_monitor`, `tavg`, `steady_tol`,
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
las paredes `wall_load_{x,y,z}` en Pa). Al final imprime el promedio
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

Split OSS (Codina): `τ_M ρ((∇u)u^k − Π[(∇u^n)u^n], (∇v)u^k)`,
`(τ_M/ρ)(∇p − Π[∇p^n], ∇q)`, `τ_C(div u − Π[div u^n], div v)`, con
`τ_M = [(2|u^k|/h)² + (4ν^k/h²)²]^{−1/2}` (sin término 2/dt por defecto:
`-kv_tau_dt true` lo agrega) y `τ_C = μ^k + ρ|u^k|h/2`. Las proyecciones
Π son L2 sobre P1 (`OrthogonalProjection.h`) y se actualizan una vez por
paso (en el iterado hacen que Picard contraiga sólo como τ/(τ+dt)).
Coeficientes (μ, τ) en P0, exactos para u en P1.

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
