# MA4703 — Resultados globales del proyecto FitzHugh–Nagumo

**Facultad de Ciencias Físicas y Matemáticas, Universidad de Chile · Primavera 2026**  
**Consolidación de resultados:** 15 de noviembre de 2026  
**Estado:** implementación numérica de las etapas A–F completada; validación de tests aprobada; optimalidad de candidatos PMP y convergencia de HJB **no certificadas**.

> **Cómo utilizar este documento.** Los enlaces a imágenes y videos son relativos a `docs/` y funcionan una vez que el archivo se copia en el repositorio. Las cifras provienen del ZIP `COPTI_Resultados_Para_Revision.zip`, generado con `scripts/recopilar_resultados.jl`, sus resúmenes JLD2, CSV, figuras y pruebas.

---

## 1. Resumen

El proyecto estudia osciladores individuales FitzHugh–Nagumo, redes de neuronas acopladas, estrategias de sincronización y desincronización, formulaciones de control óptimo, condiciones necesarias de Pontryagin (PMP) y aproximaciones de Hamilton–Jacobi–Bellman (HJB).

**Resultados centrales:**

1. Con `a=0.7, b=0.5, τ=12.5, I=1.4`, la neurona individual tiene un equilibrio único `c=(0,1.4)`, inestable, y se identificó numéricamente un ciclo periódico estable de período cercano a `34.5696`.
2. En redes de `N=12`, los acoplamientos negativos `λ=-1.5` y `λ=-0.3` produjeron convergencia numérica de estados; `λ=0.45` y `λ=1.5` mostraron dispersión importante.
3. La cota de sincronización exponencial de Lyapunov fue verificada numéricamente en los ejemplos certificados. Esta condición es **suficiente, no necesaria**.
4. En `N=4`, el control variable **factible** de sincronización alcanzó `E≤0.03` con costo `J≈14.6127`, frente a `J≈21.6627` del control constante `λ=-1.5`, para las condiciones comparadas.
5. Para desincronización con `N=4`, el candidato de corrientes balanceadas redujo `J` de `92.00` (`u=0`) a `30.2381`, con `r²(T)≈0.000636`.
6. Los residuos de PMP permanecen importantes y las aproximaciones HJB presentan cobertura insuficiente; estos resultados **no constituyen controles óptimos certificados**.

### Estado de implementación y tests

| Etapa | Contenido | Tests | Estado de implementación |
|---|---|---|
| A | EDO individual, equilibrio, estabilidad, nulclinas y retratos de fase |Implementada y probada |
| B | Hopf, ciclos, shooting, Floquet y continuación | Implementada y probada |
| C | Red completa, sincronía de estados y fase, visualizaciones | Implementada y probada |
| D | Controles suficientes, Lyapunov, corrientes diferenciadas y costos | Implementada y probada |
| E | Candidatos de control, transcripción directa y diagnóstico PMP | Implementada y probada, condiciones PMP no refinadas |
| F | HJB reducida/4D, políticas discretas y gráficos de diagnóstico | Implementada y probada, aproximación HJB sin cobertura suficiente |

> La aprobación de 445 tests confirma las comprobaciones **incluidas en el repositorio**.

## 2. Modelo y convenciones matemáticas

Modelo individual:

$$
\dot v=v-\frac{v^3}{3}-w+I,\qquad
\dot w=\frac{v+a-bw}{\tau}.
$$

Red idéntica todos-con-todos, con `z_i=(v_i,w_i)` y `\bar v=N^{-1}\sum_i v_i`:

$$
\dot v_i=v_i-\frac{v_i^3}{3}-w_i+I_0+u_i(t)+\lambda(t)(v_i-\bar v),
\qquad
\dot w_i=\frac{v_i+a-bw_i}{\tau}.
$$

**Convención de signo del enunciado:** `λ<0` es acoplamiento atractivo.

Métricas implementadas:

$$
E(x)=\frac1N\sum_{i=1}^N\lVert z_i-\bar z\rVert^2,\qquad
V(x)=\frac12\sum_i(v_i-\bar v)^2+\frac\tau2\sum_i(w_i-\bar w)^2.
$$

Para un centro fijo `c`, cuando todos los `z_i≠c`:

$$
q_i=\frac{z_i-c}{\lVert z_i-c\rVert},\qquad
r^2(x)=\left\lVert\frac1N\sum_i q_i\right\rVert^2.
$$

- `E=0` caracteriza **igualdad de estados**.
- `r=1` indica **alineamiento de las direcciones geométricas**.
- `r≈0` puede corresponder a **clústeres opuestos**, no necesariamente a ausencia de organización.
- Una caída de `E` **no garantiza** conservación de la oscilación ni alineamiento de fases, especialmente cerca de `c`.

### Problemas óptimos adoptados (Etapas E–F)

**Sincronización** — acoplamiento controlado y tiempo final libre:

$$
\min_{T>0,\,\lambda(\cdot)\in\Lambda}
J_{\rm sin}=\int_0^T\left[\rho+\alpha E(x(t))+
\beta(1-r^2(x(t)))^2+\gamma(\lambda(t)-\lambda_{\rm ref})^2\right]dt,
\quad E(x(T))\le\delta.
$$

**Desincronización** — corrientes individuales, tiempo fijo y acoplamiento fijo:

$$
\min_{u(\cdot)} J_{\rm des}=
\int_0^T\left[\alpha_d r^2(x(t))+
\gamma_d\sum_i u_i(t)^2\right]dt+\kappa r^2(x(T)),
\quad |u_i(t)|\le U_{\max},\ \sum_i u_i(t)=0.
$$

Se conservó la penalización suave **`β(1−r²)²`**, no `β(1−r)²`.

## 3. Parámetros y reproducibilidad

| Grupo de experimentos | N | Parámetros principales | Horizonte / observaciones |
|---|---:|---|---|
| A–B individual | 1 | `a=.7, b=.5, τ=12.5, I=1.4` | En B, integración y shooting de ciclos |
| C | 12 | Mismos parámetros; `λ=-1.5,-0.3,0.45,1.5` | `t∈[0,240]`, red completa; mismo estado inicial en comparaciones principales |
| D sincronización | 12 | `λ=-1.5,-0.3,-0.15,0.45`, perfil modulado y realimentación | `t∈[0,240]`, comparación de llegada a tolerancia |
| D ruptura de sincronía | 12 | `λ₀=-0.15`, corrientes balanceadas con amplitud `0,.25,.75` | `t∈[0,240]`, partida sincronizada |
| E sincronización | 4 | `δ=.03`, 10 segmentos de control `λ(t)` | Tiempo `T` variable; verificación independiente Tsit5 |
| E desincronización | 4 | `λ₀=-.15`, `|u_i|≤.8`, `Σu_i=0` | `T=28`, 10 segmentos; varias semillas iniciales |
| F sincronización | 2 | Reducción exacta especial de 2D, `δ=.03` | Iteración de Bellman sobre malla `51×41` |
| F desincronización | 2 | Cuatro estados, 5 acciones discretas `q` | `T=1.2`, malla 4D `9×7×9×7` |

Las integraciones de las etapas C–D registran `Tsit5` y tolerancias relativas/absolutas `1e-9`; B registra sus propias tolerancias `1e-9`. La comparación de costos **solo es válida entre controles del mismo experimento** (parámetros, pesos, restricciones, estado inicial y horizonte o criterio de parada comparables). En particular, los costos de `N=12`, `N=4` y `N=2` **no se comparan entre sí**.

## 4. Etapa A — Dinámica individual

Para los parámetros base:

- **Equilibrio único:** `(v*,w*)=(0,1.4)`.
- **Clasificación local:** inestable. En el equilibrio, `tr(J)=1−b/τ=0.96>0` y `det(J)=(1−b)/τ=0.04>0`.
- Las series temporales y el retrato de fase muestran comportamiento oscilatorio en el régimen analizado.

**Figuras:** [A_01: series temporales](figures/A_01_series_temporales.png), [A_02: retrato de fase](figures/A_02_retrato_fase.png), [A_03: estabilidad versus corriente](figures/A_03_estabilidad_vs_corriente.png).

## 5. Etapa B — Bifurcaciones de Hopf y ciclos límite

| Magnitud | Valor |
|---|---:|
| Corriente crítica inferior `I₋` | `0.1066694158` |
| Corriente crítica superior `I₊` | `2.6933305842` |
| Coeficiente cúbico radial `κ` | `−0.0025510204` |
| Clasificación local | **Hopf supercrítica** (`κ<0`) |
| Frecuencia angular lineal de Hopf `ω_H` | `0.28` |
| Período lineal de Hopf `2π/ω_H` | `22.4399475` |
| Período del ciclo refinado en `I=1.4` | `34.5696068` |
| Amplitud de `v` en ciclo refinado | `3.9237712` |
| Amplitud de `w` en ciclo refinado | `1.7624447` |
| Residuo de cierre por shooting | `6.23×10⁻¹⁵` |
| Multiplicadores de Floquet | aprox. `0` y `1.00000000005` |
| Clasificación de estabilidad orbital | **Estable** |

**Interpretación:** el multiplicador cercano a uno corresponde a la dirección tangente de un ciclo de un sistema autónomo, el otro muestra fuerte contracción transversal. El período `34.57` no tiene por qué coincidir con el período lineal `22.44` porque `I=1.4` no está cerca de un punto de Hopf.

**Figuras:** [B_01: equilibrios y Hopf](figures/B_01_hopf_equilibrios.png), [B_02: ciclo refinado](figures/B_02_ciclo_refinado.png), [B_03: amplitud y período](figures/B_03_amplitud_periodo.png), [B_04: continuación de ciclos](figures/B_04_continuacion_ciclos.png).

## 6. Etapa C — Redes acopladas

Experimento `N=12`, mismo estado inicial de fases distribuidas sobre un ciclo refinado, horizonte `T=240`, sin corrientes diferenciadas (`u_i=0`).

| Acoplamiento | `E(0)` | `E(240)` | `r(240)` | Observación |
|---|---:|---:|---:|---|
| `λ=-1.5` | 2.386743 | `1.19×10⁻¹⁸` | `1.000000` | Sincronización de estados numérica |
| `λ=-0.3` | 2.386743 | `1.58×10⁻¹¹` | `1.000000` | Sincronización en este experimento, fuera del criterio global suficiente |
| `λ=0.45` | 2.386743 | `4.191402` | `0.012670` | Dispersión y fuerte cancelación de fase global |
| `λ=1.5` | 2.386743 | `7.370928` | `0.007552` | Mayor dispersión en este caso |

`E(240)` muy pequeño no equivale a sincronización **exacta** en un tiempo finito; los valores están sujetos a precisión numérica. Los valores bajos de `r` tampoco demuestran inexistencia de clústeres.

**Figuras:** [C_01: acoplamientos](figures/C_01_comparacion_acoplamientos.png), [C_02: series neuronales](figures/C_02_series_temporales.png), [C_03: plano de fase](figures/C_03_retratos_fase.png), [C_04: métricas](figures/C_04_metricas.png), [C_05: grafo inicial](figures/C_05_grafo_inicial.png), [C_06: grafo final](figures/C_06_grafo_final.png).

**Videos disponibles:** `videos/C_01_red_sincronizacion.mp4` y `videos/C_02_red_sincronizacion.mp4`.

## 7. Etapa D — Controles suficientes y comparación con corrientes

### 7.1. Criterios de Lyapunov

Para `b>0`, corriente común y `λ(t)≤−1−ε`, la teoría demuestra

$$
V(t)\le V(0)e^{-2ct},\qquad c=\min\{\varepsilon,b/\tau\}.
$$

El ejemplo `λ=-1.5` y el control de corrientes `u_i=-k(v_i-\bar v)` con parámetros de la prueba verificaron la cota, el **exceso máximo de la cota informado fue 0**. Se registró `c=0.04`. La cota da un **tiempo suficiente**, no el óptimo. Para `δ=0.01`, el tiempo suficiente informado fue `100.01035`, frente a la llegada observada `20.52211` con `λ=-1.5`.


### 7.2. Sincronización para controles prefijados, `N=12`

| Estrategia | `E(240)` | Primera llegada a `δ=0.01` | Costo `J` hasta llegada | Interpretación |
|---|---:|---:|---:|---|
| `λ=-1.5` | `1.19×10⁻¹⁸` | `20.5221` | `41.0705` | Cumple el criterio suficiente de Lyapunov |
| `λ=-0.3` | `1.58×10⁻¹¹` | `105.7416` | `136.6729` | Sincroniza numéricamente |
| `λ=-0.15` | `0.101865` | No alcanzada en `[0,240]` | No definido (`+∞` en el reporte) | Cerca en fase al final, pero no alcanza la tolerancia de estados |
| `λ=0.45` | `4.191402` | No alcanzada en `[0,240]` | No definido (`+∞` en el reporte) | No sincroniza |
| `λ(t)≤-1.35`, modulado | `5.16×10⁻¹⁹` | `13.5406` | `26.6231` | Más rápido y menor costo que `λ=-1.5` en esta comparación |

`+∞` en la tabla es el **marcador de no llegada del cálculo con ese horizonte**, no la demostración de que ningún control pueda llegar a tiempo finito.

En el caso `λ=-0.15`, se obtuvo `r(240)=0.982175` pero `E(240)=0.101865`: fase casi alineada no implica todavía la tolerancia de estados elegida. Se calculó además una integral de esfuerzo de realimentación `∫Σu_i² dt≈83.11081` para el experimento correspondiente.

### 7.3. Corrientes diferenciadas, `N=12`

Con acoplamiento fijo `λ₀=-0.15`, partida sincronizada y corrientes balanceadas `Σu_i=0`:

| Amplitud de corrientes | `E(240)` | `r(240)` | `J_des` | Conclusión |
|---|---:|---:|---:|---|
| 0 | 0 | 1.000000 | `242.0000` | Conserva sincronía exacta |
| 0.25 | `0.094481` | `0.996782` | `267.0984` | Ruptura de estados leve; alto alineamiento final |
| 0.75 | `0.483335` | `0.980662` | `503.3874` | Mayor separación de estados, pero fase final todavía alineada |


**Figuras:** [D_01: certificado por acoplamiento](figures/D_01_certificado_acoplamiento.png), [D_02: certificado por realimentación](figures/D_02_certificado_feedback.png), [D_03: comparaciones](figures/D_03_comparacion_sincronizacion.png), [D_04: corrientes feedback](figures/D_04_corrientes_feedback.png), [D_05: costos de sincronización](figures/D_05_costos_J_sincronizacion.png), [D_06: ruptura de sincronía](figures/D_06_ruptura_sincronia.png), [D_07: corrientes diferenciadas](figures/D_07_corrientes_desincronizacion.png), [D_08: costos de desincronización](figures/D_08_costos_J_desincronizacion.png), [D_09: grafo con corrientes](figures/D_09_grafo_con_corrientes.png).

**Video:** `videos/D_01_corrientes_diferenciadas.mp4`.

## 8. Etapa E — Candidatos de control óptimo y Pontryagin

### 8.1. Sincronización, `N=4`, `δ=0.03`

Los dos controles se compararon para el mismo experimento y la misma condición inicial `E(0)=1.6708328`.

| Indicador | `λ=-1.5` constante | `λ(t)` candidata factible | Diferencia relativa |
|---|---:|---:|---:|
| Primera llegada a `E≤0.03` | `10.290355` | `7.258727` | `−29.46 %` en tiempo |
| Costo hasta primera llegada | `21.662688` | `14.612657` | `−32.54 %` en costo |
| Factibilidad de la tolerancia | Sí | Sí | Ambas admisibles para este criterio |

Datos complementarios del control candidato:

- Tiempo final de optimización por segmentos: `T=7.2598645`.
- `E(T)=0.02989966` verificado con Tsit5 (frente a `δ=0.03`).
- Costo discreto RK4 hasta `T`: `14.568487`; costo verificado acumulado hasta `T`: `≈14.61717`.
- Controles por tramo:
  `[-1.4665, -1.1731, -0.5276, -0.0698, 0.0911, 0.0674, -0.1057, -0.8357, -2.1109, -2.5000]`.
- Residuo máximo de la ley proyectada de Pontryagin: **`≈0.9984`**.
- Residuo informado de la condición de Hamiltoniano terminal `H(T)=0`: **`≈1.1922`**.

** Nota:** en el archivo `etapa_E_sync.jld2`, la última muestra presenta

$$
E(T)\approx 0.029900,\qquad r^2(T)\approx 0.000185\quad(r(T)\approx0.01361).
$$

Se logra **la tolerancia terminal de estados**, pero **no** el alineamiento de fases medido por este `r`. Esta situación es posible cuando los estados están próximos entre sí pero cerca del centro geométrico de fase: la dirección normalizada puede cambiar mucho aunque la distancia absoluta entre estados sea pequeña. Por tanto, **no describir este experimento como sincronización de fase alcanzada**, ni como preservación comprobada de ciclos límite. El costo incluye una penalización de fase, pero no impone `r(T)` cercano a uno como restricción terminal independiente.

### 8.2. Desincronización, `N=4`, `T=28`

| Indicador | Sin control `u=0` | Candidato de corrientes |
|---|---:|---:|
| Costo total | `92.0000` | **`30.238066`** (verificado con Tsit5) |
| Reducción relativa de costo | — | **`67.13 %`** |
| `r²(T)` | `≈1` | **`0.00063582`** |
| `E(T)` | `0` (sincronía mantenida) | `2.3768615` |
| `max |u_i|` | 0 | `0.8` |
| `max |Σ_i u_i|` | 0 | `≈1.11×10⁻¹⁶` |

**Interpretación:** el candidato cumple las restricciones de amplitud y balance, reduce el costo y deja el orden de fase global muy bajo. Sin embargo, `r²(T)≈0` no descarta grupos de fases opuestas.

**Figuras:** [E_01: sincronización candidata](figures/E_01_sincronizacion_candidata.png), [E_02: residuos PMP de sincronización](figures/E_02_residuos_PMP_sync.png), [E_03: desincronización candidata](figures/E_03_desincronizacion_candidata.png), [E_04: residuos PMP de corrientes](figures/E_04_residuos_PMP_desync.png).

## 9. Etapa F — Hamilton–Jacobi–Bellman

Las HJB completas dependen de `2N` estados. Se ensayaron dos problemas **distintos** en `N=2`; ninguno constituye una política global verificada para la red de `N=4` de la Etapa E.

### F1. Sincronización: reducción exacta especial en 2D

Con `c=(0,1.4)` y `N=2`, la dinámica de diferencias se reduce exactamente a dos variables `(d_v,d_w)`:

$$
\dot d_v=(1+\lambda)d_v-\frac{d_v^3}{3}-d_w,\qquad
\dot d_w=\frac{d_v-bd_w}{\tau},\qquad
E=d_v^2+d_w^2.
$$

| Indicador | Aproximación denominada HJB | Referencia `λ=-1.5` |
|---|---:|---:|
| Llegada a `E≤0.03` | `0.968564` | `1.777938` |
| Costo reportado `J` | `3.410716` | `3.780949` |
| Celdas con función de valor finita | **59/2091 (2.82 %)** | No aplica |
| Valor aproximado en el estado inicial `d₀` | **`+∞`** | No aplica |
| Pasos fuera de cobertura de política | **14** | No aplica |

Aunque la simulación etiquetada «HJB» llegó antes y con menor costo, **no puede afirmarse que la mejora la produjo una realimentación HJB bien resuelta a lo largo de la trayectoria**: el valor es infinito en el estado inicial y hubo pasos fuera de cobertura.

### F2. Desincronización: malla 4D, `N=2`

Se calcula una aproximación temporal en cuatro variables, `T=1.2`, malla `9×7×9×7` y acciones individuales balanceadas `u=(q,−q)` con `q∈{-0.8,-0.4,0,0.4,0.8}`.

| Indicador | Resultado |
|---|---:|
| Celdas con valor finito en el corte temporal `t=0` | **245/3969 (6.17 %)** |
| Cobertura de toda la trayectoria | **No** |
| Control aplicado fuera de cobertura | `u=0` (*fallback* declarado) |
| Control máximo resultante | `0.0` |
| `J` de trayectoria reportada | `11.6000` |
| `J` de referencia `u=0` | `11.6000` |
| `r²(T)` | `≈1` |

Este experimento **no evidenció desincronización**, principalmente porque la política no estaba disponible de manera fiable en los puntos visitados y se usó control nulo. **No es evidencia de imposibilidad de desincronizar con HJB**, ni de optimalidad de `u=0`.

Las figuras HJB se deben presentar como **diagnósticos de cobertura y limitaciones numéricas**. El gráfico F_02 corregido permite visualizar explícitamente esa falta de cobertura en vez de aparentar una solución informativa donde los valores no están definidos.

**Figuras:** [F_01: valor reducido de sincronización](figures/F_01_hjb_sync_reducida.png), [F_02: diagnóstico 4D](figures/F_02_hjb_desync_corte_4D.png), [F_03: comparaciones exploratorias](figures/F_03_comparaciones_HJB.png).

### Videos comparativos generados desde la Etapa E

- `videos/F_01_comparacion_control_variable_vs_constante.mp4`: mismo `N=4`, estado inicial y modelo, control variable factible versus `λ=-1.5`.
- `videos/F_02_comparacion_corrientes_vs_sin_control.mp4`: mismo `N=4`, estado inicial y modelo, candidato de corrientes versus `u=0`.

**No son comparaciones HJB vs PMP**: emplean las trayectorias de la Etapa E. Cuando la rama controlada llega antes a la tolerancia y se congela visualmente, ello no representa una continuación física real tras ese tiempo. Los MP4 no estaban incluidos en el paquete de revisión; sus nombres, existencia y tamaños constan en el inventario.

## 10. Resultados demostrados, calculados.

| Afirmación | Lo que encontramos |
|---|---|
| Equilibrio individual único para parámetros admisibles | Demostración en el cuaderno |
| `λ≤−1−ε` con `b>0` implica decaimiento exponencial transversal | Condición suficiente para igualdad de estados |
| Con corriente común acotada no se alcanza sincronía de estados exacta en tiempo finito desde fuera | Justifica la tolerancia `δ>0` |
| Un control de corrientes diferenciadas rompe la igualdad inicial | No garantiza fases persistentemente desalineadas |
| La tolerancia de estados implica alineamiento de fases sin hipótesis adicionales | **No en general** Ejemplo E: `E(T)≈0.03` con `r²(T)≈0.000185` |

## 11. Figuras y videos

**Selección prioritaria:**

| Recurso | Qué ilustra | Uso recomendado |
|---|---|---|
| [B_02](figures/B_02_ciclo_refinado.png) | Ciclo límite y cierre por shooting | Justificar que el modelo individual es oscilatorio |
| [C_01](figures/C_01_comparacion_acoplamientos.png) | `E(t)` y `r(t)` frente a `λ` | Ilustrar el papel del signo del acoplamiento |
| [D_01](figures/D_01_certificado_acoplamiento.png) | `V(t)` y cota exponencial | Mostrar conexión teoría–experimento |
| [D_06](figures/D_06_ruptura_sincronia.png) | Efecto de corrientes balanceadas | Diferenciar ruptura de estados y de fase |
| [E_01](figures/E_01_sincronizacion_candidata.png) | Control variable y llegada al objetivo | Comparar con estrategia constante, **aclarar el valor final de `r²`** |
| [E_03](figures/E_03_desincronizacion_candidata.png) | Corrientes, `r²(t)` y costo | Principal resultado computacional de desincronización |
| [F_02](figures/F_02_hjb_desync_corte_4D.png) | Falta de cobertura HJB | Presentar limitación del método en 4D |
| Video `F_01_comparacion_control_variable_vs_constante.mp4` | Dos redes con igual punto de partida | Exposición visual de la mejora en costo/tiempo |
| Video `F_02_comparacion_corrientes_vs_sin_control.mp4` | Dos redes, corrientes controladas o nulas | Exposición visual del efecto de control de corrientes |


## 12. Conclusiones

1. El estudio individual identifica parámetros con dinámica oscilatoria estable y permite calcular críticamente puntos de Hopf, un ciclo límite y su estabilidad de orbita.
2. El signo del acoplamiento determina un cambio cualitativo importante: valores negativos inducen atracción transversal, mientras los positivos estudiados sostienen separación entre neuronas. La cota global de Lyapunov es conservadora, se observó sincronización incluso cuando no es aplicable.
3. La desincronización por corrientes individuales es posible sin alterar el estímulo medio. Las estrategias manuales en grupos inducen separación de estados, pero no necesariamente un mejor costo ni desalineación final fuerte.
4. La optimización numérica produjo **candidatos admisibles y mejores que controles de referencia** en dos escenarios concretos. Aunque no pruebas de optimalidad global o satisfacción precisa de PMP.
5. Las métricas de sincronización son complementarias: el error `E` mide proximidad absoluta de estados y `r` mide direcciones respecto de un centro. El experimento E de sincronización enseña que **una tolerancia pequeña de estados puede coexistir con `r` muy bajo**.
6. Las formulaciones HJB están establecidas y se implementaron prototipos pero que por la dimensionalidad limitan su uso directo.