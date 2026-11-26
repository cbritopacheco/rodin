# Qué dicen los tests sobre PTAGX

Dos juegos de parámetros (`ptagx_common.py`). Los veredictos son sobre `MODEL`; los otros aparecen en cada `verdict.md` como sensibilidad. `python3 run_all.py` imprime cada veredicto apenas termina el test.

- **MODEL** (propuesto): k_a = 0.035 s⁻¹ (de E2), κ = 10 s⁻¹, ε = 0.1, monómero de fibrina, σ = 0.
- **SURFACE** (histórico): κ = 3, ε = 1, σ = 2.0×10⁻³ m calibrada en E5. Pasa E5 y rompe E2–E4; reemplazada por la trombomodulina.
- **LEGACY**: PTAGX por defecto (k_a = 0.10, κ = 3, ε = 1, σ = 0).

## Qué cambió en el modelo, y por qué

Tres cambios, dos constantes nuevas. Con `-pg_eps_x 1 -pg_kappa_x 3 -pg_fibrin_monomer false` (banco) o `-ptag_eps_x 1 -ptag_kappa_x 3 -ptag_fibrin_monomer false -ptag_k_tm 0` (aurícula) se recupera el modelo de la ronda anterior operación por operación.

### 1. El cofactor Va en el paso Xa → trombina (ε, una constante nueva)

La X de PTAGX ya era **factor Xa**, no el zimógeno: el flujo de pared j_X = Γ_FT·k_cat lleva el factor tisular y el FVIIa, y el FX del plasma (170 nM) pierde 0.8 nM en 1000 s con la dosis de Kastrup, así que llevarlo como campo sería un campo y una K_M de Michaelis de más.

Lo que faltaba es que Xa libre activa mal la protrombina: el complejo Xa·Va·fosfolípido es órdenes de magnitud más rápido, y el Va lo hace la trombina misma. κ = 3 s⁻¹ mezclaba las dos cosas en un valor "parcialmente ensamblado". Separado:

R = [k_a θ T²/(K_A+T) + κ X θ_X] P/P_ref,  θ_X = ε + (1−ε) T/(K_A+T)

con κ = 10 s⁻¹ (giro de la protrombinasa ensamblada) y ε = 0.1. La iniciación tiene que arrancar su propio cofactor, así que es un umbral; la onda que propaga, donde θ_X → 1, no cambia. Sigue siendo cuadrática en T en el origen, así que R'(0) = 0 y una traza de trombina decae.

### 2. Monómero de fibrina advectado (k_pol, de literatura)

El corte de trombina produce un monómero soluble que el flujo transporta; solo el polímero es gel:

L f = R_F − R_pol,  ∂_t F = R_pol,  R_pol = k_pol f (f + F)

con k_pol = 10³ m³ mol⁻¹ s⁻¹ (≈10⁶ M⁻¹s⁻¹, Weisel & Nagaswami 1992). G + f + F sigue conservado por el paso nodal. La vida del monómero es 1/(k_pol f): 1 s a 1 µM, 100 s a 10 nM. Por eso **sin flujo no cambia nada** (gela donde se cortó) y con flujo sí: a nM se lo lleva la corriente antes de gelar. Ese es el umbral de fibrina bajo flujo de Colace 2012.

### 3. E apaga la protección además de encender el factor tisular (k_TM, sin ajustar)

La activación inflamatoria del endotelio no sube una ganancia genérica: sube el FT y baja la trombomodulina, el EPCR, el TFPI y el eNOS a la vez. En la aurícula está medido: el marcapaseo rápido baja la expresión de TM y TFPI en el endocardio auricular y no en el ventrículo (Yamashita 2003), y el FT sube en el endotelio auricular en FA no valvular (Nakamura 2003).

Así que el mismo E hace las dos cosas: enciende la fuente, j_X E, y apaga la protección, k_TM(1−E). Dos detalles que importan:

- **(1−Ẽ) no es una capa**: vale ~1 en toda la luz. La protección es la extensión de (1−E), o sea L̃ − Ẽ con L̃ la extensión de 1, con su propia normalización de dosis: ∫(k_TM/δ)χ(L̃−Ẽ) dV = k_TM ∫(1−E) dA.
- **La contabilidad**: la trombina que atrapa la TM se inactiva por antitrombina ~20× más rápido, o sea sale como complejo serpina. El paso nodal la paga con antitrombina 1:1, con lo cual P + T − A sigue exacto y la protección queda acotada por el pool de inhibidor.

k_TM ≈ Γ_TM k_on ~ 10⁻⁶–10⁻⁵ m/s con 10²–10³ TM/µm² y un k_on casi limitado por difusión; se usa 10⁻⁵. El eNOS no entra: actúa sobre plaquetas, que PTAGX no tiene. El TFPI quedó como opción (`-ptag_tfpi_loss`, k_dX → k_dX(1−λẼ)), default 0.

**Esto reemplaza a σ.** σ era una amplificación no acotada sobre toda superficie reactiva: por eso rompía E2 (el activador sano), E4 (los parches chicos) y E3. La pérdida de TM es sustractiva, acotada, y solo existe donde hay endotelio, así que no toca ningún test de banco y no puede romper la calibración.

## Por qué la trombina necesita convección-difusión-reacción (y la fibrina polimerizada no)

La capa de trombina junto a la pared es delgada **y se la lleva el flujo**; lo segundo es lo que hay que capturar. Pe = UL/D va de 4×10⁵ en la orejuela en estasis a 2×10⁷ en el cuerpo auricular. Sobre un parche reactivo la capa de concentración es la de Lévêque, δ = (D L_p/γ)^(1/3), y el tiempo que el fluido pasa sobre el parche es t_res = L_p/(γδ):

| γ (s⁻¹) | L_p | δ | t_res | 1/k_a | ¿inicia? |
|---|---|---|---|---|---|
| 5 | 200 µm | 12 µm | 3.3 s | 29 s | no |
| 40 | 200 µm | 6 µm | 0.8 s | 29 s | no |
| 0.4 (orejuela) | 1 cm | 105 µm | 239 s | 29 s | sí |

Ese cociente t_res/(1/k_a) **es** el umbral de cizalle de Shen 2008 y **es** la hipótesis clínica (estasis → trombo). Si se quita u·∇T, el parche inicia a cualquier cizalle y la aurícula coagula en todas partes.

La crítica que sí es válida es de malla: hay que resolver δ. En E5 (hy0 = 2 µm) se resuelve. En la aurícula con h = 0.15 mm no se resolvería a cizalle alto, pero δ ∝ γ^(−1/3) y en la orejuela δ ≈ 105 µm, comparable a h. O sea: el modelo resuelve la capa donde predice trombo y no la resuelve donde no predice nada. Eso es un límite de validez, no un error, y es también por qué la capa endocárdica δ = 0.5 mm no es un truco numérico: la pared trabeculada tiene espesor.

La fibrina es el reparto en tres: fibrinógeno soluble y repuesto por el flujo (advectado), monómero soluble y de vida corta (advectado, cambio 2), polímero anclado (inmóvil). Advectar el polímero sería incorrecto: el gel está en reposo sobre una pared rígida, y es lo que lee el arrastre de Brinkman.

## Veredictos (MODEL)

| Test | Antes (κ=3, ε=1, σ=0) | Ahora (κ=10, ε=0.1) |
|---|---|---|
| E1 CAT | FALLA, 2 de 3 (pico 238 nM) | **PASA, 3 de 3**: lag 1.61 min, pico 311 nM, ETP 1595 |
| E2 Thrombodynamics | PASA, 5 de 5 | **PASA, 5 de 5**: Tlag 1.00, Vi 55.4, Vst 26.8, CS30 1160 |
| E3 DOAC (ADAGE) | FALLA, 0 de 3 | FALLA, 0 de 3 (sin cambio) |
| E4 Kastrup (axisimétrico) | FALLA, 2 de 3; umbral 100–200 µm | **FALLA, 2 de 3; umbral 50–100 µm, el medido** |
| E5 Shen (con flujo) | FALLA, 2 de 4 | FALLA, 2 de 4: no inicia a 40 s⁻¹ ✓ ni a 70/150 ✓; tampoco a 5 s⁻¹ ✗ ni a 15/700 ✗ |

Criterios: 11 de 18 antes, **12 de 18** ahora. Pero el conteo esconde lo importante: E1 pasa entero y E4 acierta el umbral de tamaño medido, que antes erraba por un factor 2.

E1 y E2 son ahora calibración (ε en E1, k_a en E2). E3, E4, E5, el diagnóstico de la onda y el caso de 1 pM de E1 son predicción.

## E3: el test que el modelo no puede pasar, y por qué

Un DOAC se modela bajando κ (anti-Xa) o subiendo A₀ (anti-IIa). El criterio es: igualar el cociente pico/valle del Tlag medido y después mirar Vi y CS30.

| caso | Tlag/ctrl | Vi/ctrl | CS30/ctrl |
|---|---|---|---|
| κ×0.9 | 1.08 | 0.97 | 1.00 |
| κ×0.75 | 1.17 | 0.90 | 0.98 |
| κ×0.6 | 1.33 | 0.82 | 0.96 |
| κ×0.5 | 1.42 | 0.74 | 0.94 |
| κ×0.25 | 2.17 | 0.41 | 0.69 |
| κ y k_a ×0.94 | 1.00 | 0.92 | 0.73 |
| **medido, rivaroxabán** | **1.67** | **0.86** | **0.73** |
| **medido, apixabán** | **1.18** | **0.95** | **0.82** |

El patrón es claro: bajando solo κ, **CS30 casi no se mueve** (0.94–1.00 contra 0.73–0.82 medido) mientras **Vi se mueve de más**. Bajando κ y k_a juntos, CS30 acierta (0.73 exacto para rivaroxabán) y Vi también (0.92 contra 0.86), pero entonces el Tlag no se mueve (1.00 contra 1.67).

O sea: el modelo puede igualar la respuesta del **coágulo** o la del **retardo**, no las dos a la vez. La razón es estructural y es la misma de la onda autónoma: el Xa solo vive a ~68 µm del activador, así que bajar κ retrasa el arranque pero no toca la onda; y la onda la mueve k_a, que un anti-Xa no toca en esta parametrización. En la realidad sí lo toca, porque la realimentación de la trombina pasa por el FXI → FIXa → FXa. Falta ese lazo.

El dabigatrán es un caso aparte: subir la antitrombina un 5–35 % mueve CS30 (0.73 → 0.39) pero no mueve el Tlag en absoluto, porque con ε = 0.1 el arranque lo fija el flujo de Xa y el umbral del cofactor.

El cofactor no cambió este test: con κ×0.5, antes daba 1.44/0.77/0.86 y ahora 1.42/0.74/0.94.

## E4: el umbral de tamaño, que es predicción pura

Nada en E4 está ajustado: la dosis j_X = Γ_FT·k_cat es la densidad medida por Kastrup. El cofactor se calibró en E1, que es 0D y no tiene ni espacio ni parches.

| parche | antes (κ=3, ε=1) | ahora (κ=10, ε=0.1) | Kastrup |
|---|---|---|---|
| 25 µm | no inicia ✓ | no inicia ✓ | no coagula |
| 50 µm | no inicia ✓ | no inicia ✓ | no coagula |
| 100 µm | no inicia ✗ | 6.0 min | coagula < 3 min |
| 200 µm | 4.6 min ✗ | 1.7 min ✓ | coagula < 3 min |
| 400 µm | 2.6 min ✓ | 1.3 min ✓ | coagula < 3 min |
| frente | 28.7 µm/min ✓ | 30.8 µm/min ✓ | 25–35 µm/min |

El umbral pasó de 100–200 µm a **50–100 µm**, que es la ventana medida. El criterio formal sigue fallando porque el parche de 100 µm tarda 6 min y se piden menos de 3, pero eso es un error de tiempo en el caso del borde, no de umbral.

## El frente sí es una onda autónoma (`source_vs_wave.py`)

Misma cinética, dosis del activador variando ocho veces:

| j_X / j_X₀ | Tlag (min) | Vi (µm/min) | Vst (µm/min) | CS30 (µm) |
|---|---|---|---|---|
| 0.5 | 1.42 | 40.9 | 30.0 | 1085 |
| 1 | 1.00 | 55.4 | 26.8 | 1160 |
| 2 | 0.67 | 62.3 | 25.1 | 1205 |
| 4 | 0.50 | 66.3 | 24.1 | 1235 |

CS30 se mueve 6 % y Vst 20 % mientras la dosis cambia 8 veces: el coágulo final no se acuerda del activador. Lo que sí depende de la dosis es el retardo y Vi, y Vi porque se mide entre 2 y 6 min *después* del retardo, cuando el frente todavía está cerca del activador.

Lo estrecho es el margen: k_a/(k_on A₀) = 0.035/0.020 = **1.76**, apenas por encima de 1, y k_a no se puede subir porque lo fijan Vst y CS30 de E2 (con k_a = 0.04 el coágulo se va a 1780 µm contra 800–1200 medidos). O sea, un solo lazo de volumen hace de motor de la onda y de limitador del tamaño final a la vez.

## E5: la dosis de factor tisular, que es lo que ese experimento restringe

Shen et al. 2008 no publica la densidad de FT de sus parches, así que lo que el experimento restringe es el producto dosis × cinética. Calibrar σ en él gastaba un parámetro de la cinética en una propiedad no medida del experimento; calibrar j_X lo gasta en esa propiedad (`calibrate_jx.py`).

| j_X (mol m⁻² s⁻¹) | × Kastrup | inicia a 5 s⁻¹ | inicia a 40 s⁻¹ | ¿ambos criterios? |
|---|---|---|---|---|
| 4.0×10⁻¹³ | 1 | no | no | no |
| 4.0×10⁻¹² | 10 | no | no | no |
| **1.2×10⁻¹¹** | **30** | **10 s** | no | **sí** |
| 4.0×10⁻¹¹ | 100 | 4 s | 10 s | no |

Hay una dosis que reproduce el umbral de cizalle completo, o sea el mecanismo está. Pero es 30 veces la única densidad medida (Kastrup), y además inicia en 10 s cuando lo medido es ~210 s: el modelo pasa de "nunca" a "de inmediato" sin el régimen lento. Eso es diagnóstico del techo estructural de abajo: el ~210 s de Shen lo fija la acumulación de complejos *unidos a la superficie*, que θ_X no representa.

## Qué se calibró y con qué, después del cambio

| Constante | Valor | De dónde |
|---|---|---|
| κ | 10 s⁻¹ | literatura: giro de la protrombinasa ensamblada (Xa·Va·fosfolípido) |
| ε | 0.1 | **ajustada en E1** (lag/pico/ETP del trombograma); robusta en [0.03, 0.3] |
| k_a | 0.035 s⁻¹ | **ajustada en E2** (sin cambio respecto de la ronda anterior) |
| k_pol | 10³ m³ mol⁻¹ s⁻¹ | literatura (Weisel & Nagaswami 1992) |
| k_TM | 10⁻⁵ m/s | Γ_TM·k_on con 10²–10³ TM/µm²; no ajustada |
| j_X | Γ_FT·k_cat | medida en E4 (Kastrup); desconocida en E2, E5 y la aurícula |

Honestidad sobre E1: antes era predicción y fallaba el pico; ahora ε se ajusta en él, así que dejó de ser predicción. Lo que sí sigue siendo informativo es que **una** constante satisface **tres** intervalos a la vez y además recupera un comportamiento cualitativo que el modelo no tenía (el pico del trombograma depende de la dosis de FT). Quedan como predicción: el caso de 1 pM de E1, E3, E4, E5 y Q1.

## Lo que deliberadamente no se agregó

- **FVII/FVIIa y FT como campos.** El FVII del plasma está a 10 nM y el FVIIa a ~0.1 nM, los dos prácticamente constantes; el complejo FT·FVIIa es una especie *de superficie*, confinada a la pared, con ensamblaje rápido (Kd ~ 1 nM). Su densidad cuasi-estacionaria por k_cat **es** j_X. Llevarlos agregaría dos campos de volumen, una especie de superficie y tres constantes para no cambiar nada.
- **El zimógeno FX.** 170 nM que pierden 0.8 nM por corrida: un campo y una K_M de más.
- **eNOS/NO.** Actúa sobre plaquetas, que PTAGX no tiene por construcción.
- **φ, el lazo intrínseco con retardo** (FXI → FIXa → FXa como fracción retardada de la realimentación, `-pg_phi_i`). Implementado y probado en E2: φ = 0.3 da Tlag 0.83, Vi 37.1, Vst 23.5, CS30 905 (Vi se sale por 2 %) y φ = 0.7 saca a Vst y CS30. Además va en la dirección equivocada para E3: al enrutar parte de la *propagación* por Xa, la vuelve **más** sensible al anti-Xa, cuando los datos piden menos. Queda en el código con default 0 y documentado como descartado.
- **Complejos enzimáticos de superficie con su propio tiempo de residencia.** Es el techo estructural que queda (ver abajo).

## El techo que queda

θ_X es un sustituto *de volumen* de un hecho *de superficie*: la protrombinasa se ensambla sobre la membrana de la lesión y se queda ahí mientras el flujo barre el plasma. Un modelo con especies de superficie (Xa unido, Va unido, con sus tasas de adsorción y desorción) sería el paso siguiente y es lo que haría falta para la iniciación bajo flujo de E5 con la dosis medida. Costo: dos campos de superficie y cuatro constantes. No se hizo aquí a propósito.
