# ArterialLesionAxi — fluidos CY y sPTT generalizado, y jacobiano corregido (9 oct 2026)

## Qué cambió en el driver

1. **`-al_fluid sptt | cy | gsptt`** (default `sptt`, comportamiento anterior intacto).
   - `cy`: Carreau–Yasuda generalizado-newtoniano. Problema (u, p) aparte, sin ψ, un solve por paso.
     Viscosidad evaluada en uⁿ (retrasada, primer orden como BDF1), incluida la deformación de giro u_r/r.
   - `gsptt`: un modo sPTT + solvente Carreau–Yasuda con η_s(γ̇) = η_CY(γ̇) − η_p/f(γ̇).
     La viscosidad estacionaria del modelo es **exactamente** η_CY: `cy` es su gemelo sin memoria por construcción.
     Exige λ₀ = λ. El código comprueba al arrancar que η_s(γ̇) > 0 en 10⁻³–10⁵ s⁻¹.
   - Parámetros CY: `-al_cy_eta_zero 0.056 -al_cy_eta_inf 0.00345 -al_cy_lambda 1.902 -al_cy_n 0.22 -al_cy_a 1.25` (Cho y Kensey 1991, defaults).
2. **`-al_eta_ref <Pa s>`**: viscosidad que define Re y Wo. Con 0 (default) se usa η_s + η_p (sptt) o η_∞ (cy, gsptt).
   Con `-al_eta_ref 0.00345` los cuatro fluidos corren al **mismo flujo físico** sin convertir Re ni Wo a mano.
3. **Jacobiano de Newton corregido** en el término de estabilización S1 de div σ: faltaba la derivada de Dexp respecto a ψ
   (D²exp[∂ψᵏ, δψ]). Con ψ grande (sangre real) Newton convergía linealmente (8–10 iteraciones); ahora es cuadrático (3–5).
   La solución convergida no cambia; solo el número de iteraciones. Los drivers planar y 3D tienen la misma omisión.
4. **`-al_inlet_developed_conformation 1`**: impone en la entrada la conformación de corte simple estacionario del perfil
   medio de Poiseuille (ψ = log c_ss), en vez de ψ = 0. Con λ Ū/D ≫ 1 la entrada relajada necesita ~λŪ de tubo para
   desarrollarse (300 D para el fluido R a Re 300). Los scripts lo activan para R; para D (λŪ/D ≈ 2) da igual.
5. `-al_trace_newton 1`: imprime la historia de Newton en todos los pasos (diagnóstico).
6. `PointwiseScalar` acepta un orden de cuadratura: con viscosidad constante se usa orden 0 → el ensamblaje `sptt` es
   **bit a bit el anterior** (regresión: CSV idéntico a 6 cifras en 14 pasos, difiere solo el último dígito de |dψ|).

## Verificación (malla H axisimétrica rf 0.25, 3,7k vértices, Re 50 con η_ref 3,45 mPa·s, Wo 2, flujo estacionario)

| Prueba | Resultado |
|---|---|
| Poiseuille CY vs. solución semianalítica (ODE radial integrada) | Δp 53,55 Pa vs 53,76 (−0,4 %); caudal 1,4 % bajo por los 5 elementos radiales de la malla de prueba |
| gsptt con entrada relajada (ψ = 0) | Δp 52,07 Pa (−2,8 %): el polímero del núcleo no se desarrolla en 35 D (λŪ = 50 D); máx σ 5,3 Pa por el sobreimpulso de arranque del sPTT en la pared |
| gsptt con entrada desarrollada | Δp 53,76 Pa = CY dentro de 0,4 % (propiedad gemela); máx σ 2,0–2,3 Pa = σ_xx analítico en la pared (2,09) |
| Newton gsptt sin corrección del jacobiano | 8–10 iteraciones, convergencia lineal (razón ≈ 0,25) |
| Newton gsptt con corrección | 3–5 iteraciones, cuadrática (1,6 → 0,28 → 0,07 → 4,5e-4 → 1,3e-6) |
| Verificación externa del mapa 3×3 (Python) vs. código | error 1e-10 (sin cambios) |
| Piloto pulsátil S75, Re 100, Wo 2,5, fluido R (binario anterior) | estable hasta el pico del primer ciclo (σ 9 Pa, 7 iteraciones); con el jacobiano nuevo, ver `pilot_new` |

## Cambios del 9 oct, mediodía

7. **Entrada desarrollada para `cy` y `gsptt`**: el modo medio del perfil de entrada ya no es la parábola newtoniana sino el
   perfil estacionario desarrollado de Carreau–Yasuda a Ū (EDO radial, G por bisección); los modos oscilatorios siguen siendo
   Womersley con η∞ (linealización en torno al flujo medio). La conformación desarrollada de la entrada (`gsptt`) usa la
   tasa de corte de ese mismo perfil. Verificado contra Python: τ_w 0,3840 Pa, γ̇_w 86,7 s⁻¹, u_max/Ū 1,835 (idénticos).
8. **Índices de pared sobre la lesión**: el CSV agrega `maxTAWSSLesion,maxOSILesion` (máximo sobre el tag 5), al final de la
   fila para no romper scripts existentes. La línea `[cycle ...]` del log también los imprime. Para el paper usar estos,
   no `maxTAWSS` global (contaminado por la entrada en H y aneurismas).

Resultado en el tubo recto (Re 50, η_ref 3,45 mPa·s, TAWSS analítico 0,384 Pa):

| | CY antes | CY ahora | R antes | R ahora |
|---|---|---|---|---|
| TAWSS a x/D = −9,5 | 0,379 | 0,390 | 0,420 | 0,427 (instantáneo bajando 0,447 → 0,427 entre ciclos) |
| máx. TAWSS global | 0,3995 | 0,399 | 0,465 | 0,473 |
| TAWSS máx. en la lesión (nuevo) | — | 0,3827 (−0,3 %) | — | 0,3841 (0,0 %) |

En R el pico cerca de la entrada **no** venía del perfil de entrada (el diagnóstico de la mañana era incorrecto para R): persiste
con velocidad y conformación desarrolladas y consistentes, y decae lentamente entre ciclos. Lo más probable es un transitorio
de arranque local (interior en ψ = 0 al inicio, relajación con λ = 4,8 s cerca de la pared donde la advección es lenta).
Está confinado a x/D < −8 (8 D aguas arriba de la lesión) y queda excluido de los índices de la lesión.

## Scripts

- `run_axi_re.sh <Re>`: fluidos `FLUIDS="D R"` por defecto, `WI_LIST="2"`, Wo por Re (100→2,5; 300→4; 600→5,5).
  Carpetas `$OUT/Re<Re>_Wo<Wo>/<geo>/{N,Wi<wi>,CY,R}`. Reutiliza `$OUT/Re<Re>` si Wo = 4 (barrido anterior).
  R corre con `CYCLES_R=8`, `-al_conformation_its 10`, De 5,635, entrada desarrollada.
- `run_axi_extra.sh <Re>`: el mismo script con `FLUIDS="N CY"` (baratos), para otra terminal.
- Ambos escriben `summary.csv` con fluido, Wi, pasos, ciclos, estado y minutos. Reanudables.

## Pendiente / riesgos

- `cy` y `gsptt` usan η_∞ en el perfil de Womersley de entrada (no hay perfil desarrollado cerrado para CY). Entrada a 10 D.
- Para R, la periodicidad entre ciclos hay que mirarla en el CSV (λ ≈ 6 latidos); 8 ciclos es una estimación.
- Los drivers 2D planar y 3D no tienen la corrección del jacobiano ni los fluidos nuevos.
- Nada está commiteado en git.
