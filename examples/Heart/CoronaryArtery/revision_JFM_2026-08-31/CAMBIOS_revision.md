# Registro de cambios — revisión del manuscrito JFM (31-08-2026)

Archivo base: el `.tex` entregado por Oscar (versión "reoriented" con macros
editoriales). Archivo revisado: `JFM_Rmu_coronary_rev.tex`.
Política aplicada: cambios mínimos; sólo se promovió a texto definitivo lo que
pudo verificarse contra el código (`CoupledLV0DCoronary3D.h/.cpp`) o por
verificación numérica independiente. Lo no verificable quedó como comentario
rojo, ahora con una nota azul de estado.

## 1. Comentarios rojos RESUELTOS (eliminados; texto promovido)

1. **Referencia para la reducción sistólica de flujo coronario** (Intro).
   Se citó `\citep{Downey1975,Permutt1963}` en la frase sobre compresión
   intramiocárdica: ambas entradas existen en
   `examples/Heart/CoronaryArtery/JFM_nuevas_referencias.bib` y son la fuente
   ya usada para el modelo waterfall. Comentario eliminado.

2. **Verificación de las expresiones analíticas Cross / Carreau–Yasuda /
   Quemada** (§3 y Apéndices A y B). Se ejecutó una verificación numérica
   independiente (script `verify_closures.py`, incluido en esta carpeta):
   - Cross (ec. I_Cr con 2F1): discrepancia relativa máxima **3×10⁻¹³** frente
     a cuadratura adaptativa de la integral reológica, 3 estados × 10 valores
     de τ_w en 10⁻⁴–10² Pa.
   - Carreau–Yasuda (ec. I_CY, incl. términos J_{b−1}, J_{2b−1}, J_{3b−1}):
     máx. **2×10⁻¹²**.
   - Límite newtoniano (μ₀→μ∞): exacto a precisión de máquina en ambos.
   - Quemada F(α,q) con P₁–P₈ transcritos: acuerdo **<2×10⁻⁴** (limitado por
     la derivada numérica de la cuadratura de referencia), 3 estados.
   - Tangente exacta (ec. exact-derivative) vs diferencias finitas centradas:
     acuerdo **<10⁻¹⁰**.
   Los tres comentarios de verificación se sustituyeron por párrafos que
   reportan estos números. En el Apéndice B queda un comentario corto: falta
   sólo el cotejo visual final contra Popel & Enden (1993).

3. **Interpretación N-vasos en paralelo** (§5.1, §5.3). Confirmado en el
   código: la calibración calcula N_a = Q_i/(π r_a² v_a) por outlet (split de
   Murray), y R_a = ΔP_a/Q_i coincide algebraicamente con la resistencia de
   Poiseuille de N_a tubos en paralelo. Se añadió un párrafo al final de §5.1
   que lo dice, y se eliminaron los dos comentarios (§5.1 y el de la tangente
   en §5.3, cuya parte numérica quedó verificada, ver punto 2).

4. **Tabla de error de calibración (tabla 2)**. Recalculada íntegramente de
   forma independiente con los parámetros CY de la tabla 1: todos los valores
   e(q), μ_eff(q₀) y e_∞ se reproducen con la precisión publicada. Se añadió
   la frase que lo dice; el comentario quedó reducido a "regenerar sólo si
   cambia el ajuste".

5. **Topología Q_i del elemento proximal** (§6.2). Confirmado en el código:
   la caída proximal se ensambla implícitamente sobre la frontera de salida
   (R_a Φ_a A (u·n)(v·n)), por lo que el flujo del elemento proximal ES el
   flujo 3D Q_i; no falta ningún estado. Se añadió la frase y se eliminó el
   comentario.

6. **α, β, p_w y forma implementada de p_out** (§6.2.1). Del código:
   - α = 0.7 (`intramyocardialFraction`), promedio ponderado transmural;
     texto añadido con las referencias waterfall.
   - El código implementa p_drain = max(p_im, P_RA), es decir **β = 1**. Se
     eliminó el parámetro β y la ec. del waterfall quedó
     p_w = max(p_im, P_RA) (una ecuación y un parámetro menos — esto además
     reduce las ecuaciones de esa sección, como pedía Oscar).
   - Gate de Starling: cuando p_im > P_RA la válvula venosa se cierra al
     llegar p_tm→0; en el código está regularizada con una conductancia
     residual 10⁻³ y con la cota p_tm ≥ 0 impuesta. Texto añadido.
   - La ec. p_out se reescribió en la forma implementada: caída proximal por
     resistencia dependiente del flujo R_{μ,p}(Q_i) (parametrización en
     flujo de la misma curva monótona), NO por el mapa inverso. Comentario
     eliminado.

7. **Protocolo de simulación** (§6). Confirmado del código y del PDF:
   CY específico del estado en dominio 3D y outlets; en el caso constante
   sólo cambia el cierre de salida; sin barrido de P_RA. Se añadió la frase
   "No sweep of the right-atrial pressure is included..." y se eliminó el
   comentario.

8. **Parámetros ventriculares y sistémicos** (§6.1). Se añadió la tabla
   compacta pedida (nueva tabla `tab:lv-parameters`) con los valores del
   código: R₀=24.5 mm, d₀=14 mm, E_s=3×10⁶ Pa, σ₀=1.2×10⁵ Pa, k₀=10⁵ Pa,
   μ=η=70 Pa·s, R_p=5×10⁷, C_p=5×10⁻⁹, R_d=10⁸, C_d=5×10⁻¹⁰, p_sv=1000 Pa,
   T=0.85 s, presión auricular 500–1250 Pa. Se declaró: 3 ciclos calculados,
   2 descartados, promedios sobre el tercero. Comentario eliminado.

9. **Discretización numérica** (§6.3). Promovidos del código:
   ε_γ = 10⁻³ s⁻¹; estabilización de backflow (multiplicador 1) agrupada con
   la corrección de Temam (como pedía el comentario del PDF); impedancia
   normal de entrada 5×10² Pa·s/m; penalización tangencial 10³;
   regularización de presión 10⁻¹² con su justificación (fija el modo de
   presión para la factorización directa; Taylor–Hood ya es inf-sup);
   Δt = 10⁻³ s con adaptividad (mitades, ≤8 niveles); solver directo MUMPS
   LU (PETSc); Newton escalar del outlet: tol 10⁻⁹ Pa, ≤50 iter; tabla WRMS
   de 241 nodos log en τ_w∈[10⁻⁶,10⁴] Pa, error de interpolación máx 0.09%.
   Se afirmó que el set CY 3D es idéntico entre ambos cierres a estado fijo.
   Quedan DOS comentarios cortos: independencia de malla e independencia
   temporal (siguen pendientes de verdad — no se pueden resolver sin correr).

10. **Constantes de calibración de outlets** (§6.4). Añadidas del código:
    Q_tot = 1.5×10⁻⁶ m³/s repartido por ley de Murray (Q_i ∝ r_i³),
    f_v = 0.13, C_tot = 4×10⁻¹⁰ m³/Pa, μ_ref = 3.5 mPa·s.

11. **Figura 5 (esquema)**: sustituido el placeholder por la nueva figura
    `coronary_setup_rcrmu.pdf` (ver §3 abajo) con caption nuevo que describe
    las condiciones de borde. Comentario de redibujo eliminado.

## 2. Comentarios rojos ACTUALIZADOS (siguen abiertos, ahora con datos)

- **Tabla 1 / Quemada**: evidencia fuerte de que las filas "k₀, k∞" son en
  realidad μ₀, μ∞ en mPa·s (coinciden con los plateaus de Cross/CY del mismo
  estado). Invirtiendo la ec. (quemada_limits) con los μ_F tabulados:
  (k₀,k∞) = (4.05,1.91), (3.89,1.63), (3.62,1.15) — todos admisibles y
  consistentes con las correlaciones tipo Cokelet del solver (k₀≈4.2,
  k∞≈1.8 a φ=0.42). El comentario ahora pide relabel o reporte de los k
  ajustados. Unidades de m (power law): sigue abierto.
- **P_RA**: el manuscrito y el PDF dicen 1000 Pa, pero el default actual de
  `CoupledLV0DCoronary3D.h` es **1800 Pa**. El comentario ahora pide
  reconciliar cuál valor generó los resultados publicados.
- **Ajuste poblacional (expectation fit) y Krieger–Dougherty**: los scripts
  de ajuste NO están en el repositorio rodin (sólo `nonnewtonian_comparison.py`
  y `verify_outlet_0d.py`, que no ajustan). Se añadió nota de estado; sigue
  abierto. Nota: para Quemada el solver usa correlaciones exponenciales tipo
  Cokelet, consistente con la forma exponencial de (phi_scaling).
- **Referencias de radios/velocidades microvasculares**: se confirmó que los
  valores del texto coinciden con la implementación (r_a=25 μm, v_a=5 mm/s,
  r_v=30 μm, v_v=3 mm/s); falta la cita de microscopía intravital.
- **Gap de literatura (modelos reducidos no newtonianos)**: se sugirió en
  azul Ghigo, Lagrée & Fullana (JNNFM 253, 2018) + Sochi (2015) como
  antecedentes más cercanos; confirmar y añadir a jfm.bib.

## 3. Figura nueva

`coronary_setup_rcrmu.pdf` (fuente TikZ: `coronary_setup_rcrmu.tex` +
`mesh_cropped.pdf`, recorte del `coronaria_esquema` de Gmsh):
- fondo blanco (el comentario "transparent" del PDF queda atendido);
- malla coronaria con los 6 outlets marcados Γ_μ^(i)–R_μRCR;
- condiciones de borde: tracción −p_ar(t)n en Γ_in, no-slip en Γ_w;
- bloque LV-0D (activación → LV Caruel 2014 → Windkessel sistémico) con
  p_ar(t) hacia el inlet y p_LV(t) hacia p_im = α p_LV (línea discontinua);
- inset del outlet terminal: R_{μ,p}(Q_i) y R_{μ,d} como resistencias
  VARIABLES (símbolo con flecha), compliance C^(i) referida a p_im(t),
  drenaje waterfall p_w = max(p_im, P_RA) hacia P_RA, con las ecuaciones
  del modelo de un estado.
Reproduce exactamente (6.7)–(6.10) en su forma final (sin el capacitor+
fuente superado), como pedía el comentario.

## 4. Comentarios "sencillos" del PDF incorporados

- p1: "In this work, we derive…" y "Results show that…" en el abstract
  ("is tested" ya estaba).
- p4: "shear-thinning behaviour" ya estaba en el texto base.
- p5: "In this work, the analysis is restricted to purely viscous…" (§2.2).
- p8: power law introducido como "By writing … we infer the relation" (§3.1).
- p10: fibrinógeno como mediador de rouleaux (§2.2); φ definida como
  fracción de volumen de RBC en plasma (§2.2.4, Quemada).
- p12: asimetría del intervalo por no linealidad en φ — ya estaba.
- p14 ("transición abrupta / miedo a repetición"): la apertura de §5 se
  reescribió como transición ("The closure of §3 can now be inserted…")
  en lugar de re-introducir los modelos Windkessel; se eliminó además una
  frase redundante al final del párrafo 5 de la introducción ("displacement
  along the constitutive curve changes the effective hydraulic resistance…",
  que repetía la idea que reaparece en §5.2 y §7).
- p24: backflow stabilization mencionada y agrupada con el Temam trick en la
  misma frase de estabilización (§6.3).
- p21 ("transparent"): resuelto con la nueva figura en fondo blanco.
- p31 ("3D fluid"): no se pudo determinar sin ambigüedad a qué frase
  apuntaba la nota; no se cambió nada (avisar si se refería a otra cosa).

## 5. Otros

- Nota de bibliografía del preámbulo actualizada: Downey1975 y Permutt1963
  están en `JFM_nuevas_referencias.bib` y resuelven en el build PROPOSAL;
  hay que fusionarlos en `jfm.bib` de este archivo. Los otros tres keys
  también resuelven en el PROPOSAL, así que las entradas existen.
- No se tocó: tabla 1 (valores), resultados numéricos (§7), conclusiones
  (salvo nada), apéndices (salvo los párrafos de verificación).

## 6. Verificación numérica (resumen del script)

Ver `verify_closures.py` y `verify_closures_output.txt`. Además de lo del
punto 1.2: los γ_tr de transición recalculados dan 45/24/104 s⁻¹ (ε=0.2) y
108/55/253 s⁻¹ (ε=0.1) frente a 47/25/108 y 110/56/257 del texto —
consistentes dentro de la tolerancia de implementación (~3–5%), sin cambio
requerido.

---

## 7. Parche de código para los escenarios patológicos (01-09-2026)

Respaldos: `*.prepatch_pra` (primer parche) y `*.prepatch_const` (segundo).
Diff acumulado: `patch_operating_pra.diff`.

**7.1 `Config::operatingRightAtrialPressure`** (default 0 = comportamiento
actual). Separa la presión de drenaje *en operación* de la que usa la
calibración en reposo: un escenario de hipertensión venosa desplaza el punto
de operación de un lecho cuya (R_a, R_v, C) sigue siendo la del baseline
sano. Consumida en `updateOutlet0D` por `pDrain` y por el gate del waterfall.

**7.2 `Config::constantOutletResistance`** (default false). Congela el cierre
0D en la meseta de alto corte de la ley activa — mu_inf para Carreau-Yasuda,
mu_F(1-k_inf phi/2)^-2 para Quemada — de modo que la tabla WRMS se construye
con mu(gamma) == mu_inf constante. Con mu constante la integral reológica es
I = tau_w^4/(4 mu) *exactamente*, así que mu_ap == mu_inf en cada nodo y cada
salida se vuelve el resistor lineal R_inf = 8 mu_inf L/(N pi R^4). Clave: NO
toca `Config::viscosity`, que sigue gobernando el campo Carreau-Yasuda 3D —
por eso el par de corridas aísla el cierre reducido, que es la comparación de
§7.5. (Poner `viscosity.mu0 = viscosity.muInf` habría vuelto newtoniano
también el dominio resuelto: es el error que este flag evita.)

**7.3 Opciones del driver** (`CoronaryArtery.cpp`):
`-coronary_operating_pra <Pa>`, `-coronary_alpha_im <0..1>`,
`-coronary_compliance_total <m^3/Pa>`, `-coronary_constant_outlet <0|1>`,
`-coronary_output_prefix <dir>` (redirige CSV *y* XDMF/HDF5, así que cada
escenario queda aislado y no pisa `hyp2/`).

**7.4 Scripts** (40 hilos OpenMP por caso, corridas secuenciales, sin
requerir mpirun, se ejecutan desde `build/`; flags PETSc explícitos en
`PETSC_OPTS`, idénticos a los que se usan interactivamente):
- `run_scenarios.sh` — familia R_mu: pra2200, alpha085, combined. El caso
  baseline NO se vuelve a correr: se toma la solución de referencia ya
  existente vía `BASELINE_CSV` (default `hyp2/CoronaryArtery.csv`), que además
  fija el Q0 con el que se adimensionaliza gamma en todos los escenarios. Para
  recalcularlo, descomentar `run baseline`.
- `run_scenarios_constant.sh` — gemelos con `-coronary_constant_outlet`,
  salidas en `*_const/`, mismo Q0 de referencia, baseline igualmente omitido.
- `retrograde_metrics.py` — por ciclo: fracción de volumen retrógrado,
  duración/pico/episodios de la fase de reversión, ratio sístole/diástole,
  gamma_d (mín, media, % del ciclo bajo umbral), p_tm mínimo, extremos de Phi.
- `compare_closures.py` — comparación pareada R_mu vs constante a estado
  fijo: d<|Q_in|>, d<mu_eff> por elemento, E_Q(t) medio y extremo con su fase,
  y max mu_eff/mu_inf. Alinea las dos corridas por fase del ciclo e interpola,
  así que una malla temporal distinta (adaptividad) no sesga la comparación.
  Produce directamente las cifras de §7.5.
