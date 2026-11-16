# MA4703 — Resultados globales del proyecto FitzHugh–Nagumo

**Facultad de Ciencias Físicas y Matemáticas, Universidad de Chile · Primavera 2026**  
**Consolidación de resultados:** 15 de noviembre de 2026  
**Estado:** implementación numérica de las etapas A–F completada; validación de tests aprobada; optimalidad de candidatos PMP y convergencia de HJB **no certificadas**.

> **Cómo utilizar este documento.** Este archivo está pensado para guardarse como `docs/RESULTADOS_GLOBALES.md` dentro de `Proyecto-COPTI`. Los enlaces a imágenes y videos son relativos a `docs/`; funcionan una vez que el archivo se copia en el repositorio. Las cifras provienen del ZIP `COPTI_Resultados_Para_Revision.zip`, generado con `scripts/recopilar_resultados.jl`, sus resúmenes JLD2, CSV, figuras y pruebas. Los videos se inventariaron, pero **no fueron incluidos ni inspeccionados directamente** en ese ZIP. Se mantienen separados los resultados demostrados, los experimentos numéricos y las cuestiones pendientes.

---

## 1. Resumen ejecutivo

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
|---|---|---:|---|
| A | EDO individual, equilibrio, estabilidad, nulclinas y retratos de fase | 61/61 | Implementada y probada |
| B | Hopf, ciclos, shooting, Floquet y continuación | 38/38 | Implementada y probada |
| C | Red completa, sincronía de estados y fase, visualizaciones | 89/89 | Implementada y probada |
| D | Controles suficientes, Lyapunov, corrientes diferenciadas y costos | 161/161 | Implementada y probada |
| E | Candidatos de control, transcripción directa y diagnóstico PMP | 54/54 | Implementada y probada; condiciones PMP no refinadas |
| F | HJB reducida/4D, políticas discretas y gráficos de diagnóstico | 42/42 | Implementada y probada; aproximación HJB sin cobertura suficiente |
| **Total** | **Seis etapas** | **445/445** | **Todas las pruebas programadas aprobadas** |

> La aprobación de 445 tests confirma las comprobaciones **incluidas en el repositorio**, no la optimalidad de los controles ni la convergencia de los esquemas HJB. El recopilador ejecutó los tests A–F y pudo leer los 18 archivos JLD2 inventariados sin errores de lectura.

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

**Convención de signo del enunciado:** `λ<0` es acoplamiento atractivo; no confundir con formulaciones que emplean el Laplaciano con signo contrario.

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
- Las series temporales y el retrato de fase muestran comportamiento oscilatorio en el régimen analizado; el estudio detallado de ciclos corresponde a B.

**Figuras:** [A_01: series temporales](figures/A_01_series_temporales.png), [A_02: retrato de fase](figures/A_02_retrato_fase.png), [A_03: estabilidad versus corriente](figures/A_03_estabilidad_vs_corriente.png).

**Lectura científica:** la existencia del equilibrio y su clasificación provienen de la teoría del modelo; las curvas muestran una realización numérica concreta.

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

**Interpretación:** el multiplicador cercano a uno corresponde a la dirección tangente de un ciclo de un sistema autónomo; el otro muestra fuerte contracción transversal. El período `34.57` no tiene por qué coincidir con el período lineal `22.44` porque `I=1.4` no está cerca de un punto de Hopf.

Se continuaron **cinco** puntos entre `I=1.4` e `I=1.2`, todos con estado `converged`. Esto es continuación local con corrector por shooting, **no** un seguimiento global garantizado de todas las ramas periódicas. En un barrido se detectaron 13/17 candidatos de órbita; no constituyen automáticamente 13 ciclos distintos.

**Figuras:** [B_01: equilibrios y Hopf](figures/B_01_hopf_equilibrios.png), [B_02: ciclo refinado](figures/B_02_ciclo_refinado.png), [B_03: amplitud y período](figures/B_03_amplitud_periodo.png), [B_04: continuación de ciclos](figures/B_04_continuacion_ciclos.png).

## 6. Etapa C — Redes acopladas

Experimento `N=12`, mismo estado inicial de fases distribuidas sobre un ciclo refinado, horizonte `T=240`, sin corrientes diferenciadas (`u_i=0`).

| Acoplamiento | `E(0)` | `E(240)` | `r(240)` | Observación |
|---|---:|---:|---:|---|
| `λ=-1.5` | 2.386743 | `1.19×10⁻¹⁸` | `1.000000` | Sincronización de estados numérica |
| `λ=-0.3` | 2.386743 | `1.58×10⁻¹¹` | `1.000000` | Sincronización en este experimento, fuera del criterio global suficiente |
| `λ=0.45` | 2.386743 | `4.191402` | `0.012670` | Dispersión y fuerte cancelación de fase global |
| `λ=1.5` | 2.386743 | `7.370928` | `0.007552` | Mayor dispersión en este caso |

**Precaución:** `E(240)` muy pequeño no equivale a sincronización **exacta** en un tiempo finito; los valores están sujetos a precisión numérica. Los valores bajos de `r` tampoco demuestran inexistencia de clústeres.

**Figuras:** [C_01: acoplamientos](figures/C_01_comparacion_acoplamientos.png), [C_02: series neuronales](figures/C_02_series_temporales.png), [C_03: plano de fase](figures/C_03_retratos_fase.png), [C_04: métricas](figures/C_04_metricas.png), [C_05: grafo inicial](figures/C_05_grafo_inicial.png), [C_06: grafo final](figures/C_06_grafo_final.png).

**Videos disponibles en el repositorio:** `videos/C_01_red_sincronizacion.mp4` y `videos/C_02_red_sincronizacion.mp4`. El segundo archivo existe en el inventario; su contenido específico debe comprobarse viendo el MP4, no inferirse del nombre.

## 7. Etapa D — Controles suficientes y comparación con corrientes

### 7.1. Criterios de Lyapunov

Para `b>0`, corriente común y `λ(t)≤−1−ε`, la teoría demuestra

$$
V(t)\le V(0)e^{-2ct},\qquad c=\min\{\varepsilon,b/\tau\}.
$$

El ejemplo `λ=-1.5` y el control de corrientes `u_i=-k(v_i-\bar v)` con parámetros de la prueba verificaron la cota; el **exceso máximo de la cota informado fue 0**. Se registró `c=0.04`. La cota da un **tiempo suficiente**, no el óptimo. Para `δ=0.01`, el tiempo suficiente informado fue `100.01035`, frente a la llegada observada `20.52211` con `λ=-1.5`.

Bajo campo medio, `u_i=-k(v_i−\bar v)` y acoplamiento basal `λ₀` tienen dinámica equivalente a acoplamiento efectivo `λ₀−k`, mientras no se agreguen otras corrientes no equivalentes. Se comprobó la consistencia de esta relación con pruebas.

### 7.2. Sincronización para controles prefijados, `N=12`

| Estrategia | `E(240)` | Primera llegada a `δ=0.01` | Costo `J` hasta llegada | Interpretación |
|---|---:|---:|---:|---|
| `λ=-1.5` | `1.19×10⁻¹⁸` | `20.5221` | `41.0705` | Cumple el criterio suficiente de Lyapunov |
| `λ=-0.3` | `1.58×10⁻¹¹` | `105.7416` | `136.6729` | Sincroniza numéricamente; no cumple el criterio global usado |
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

Las corrientes con perfil *split* imponen grupos con estímulos diferentes. La aparición visual de dos grupos **no equivale por sí sola** a desincronización de fase sostenida. Tampoco hubo mejora del costo respecto de `u=0` en estos tres ensayos manuales.

**Figuras:** [D_01: certificado por acoplamiento](figures/D_01_certificado_acoplamiento.png), [D_02: certificado por realimentación](figures/D_02_certificado_feedback.png), [D_03: comparaciones](figures/D_03_comparacion_sincronizacion.png), [D_04: corrientes feedback](figures/D_04_corrientes_feedback.png), [D_05: costos de sincronización](figures/D_05_costos_J_sincronizacion.png), [D_06: ruptura de sincronía](figures/D_06_ruptura_sincronia.png), [D_07: corrientes diferenciadas](figures/D_07_corrientes_desincronizacion.png), [D_08: costos de desincronización](figures/D_08_costos_J_desincronizacion.png), [D_09: grafo con corrientes](figures/D_09_grafo_con_corrientes.png).

**Video:** `videos/D_01_corrientes_diferenciadas.mp4`.

## 8. Etapa E — Candidatos de control óptimo y Pontryagin

### 8.1. Sincronización, `N=4`, `δ=0.03`

Los dos controles se compararon para el mismo experimento y la misma condición inicial `E(0)=1.6708328`. La trayectoria de control variable utiliza 10 segmentos.

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

**Advertencia central de interpretación de fase:** en el archivo `etapa_E_sync.jld2`, la última muestra presenta

$$
E(T)\approx 0.029900,\qquad r^2(T)\approx 0.000185\quad(r(T)\approx0.01361).
$$

Se logra **la tolerancia terminal de estados**, pero **no** el alineamiento de fases medido por este `r`. Esta situación es posible cuando los estados están próximos entre sí pero cerca del centro geométrico de fase: la dirección normalizada puede cambiar mucho aunque la distancia absoluta entre estados sea pequeña. Por tanto, **no describir este experimento como sincronización de fase alcanzada**, ni como preservación comprobada de ciclos límite. El costo incluye una penalización de fase, pero no impone `r(T)` cercano a uno como restricción terminal independiente.

Se trata de un **control factible y comparativamente mejor** que el control constante ensayado, pero sus residuos PMP no permiten calificarlo como extremal normal verificado ni como óptimo global. El ligero margen `E(T)<δ` tampoco certifica actividad exacta de la frontera para un óptimo.

### 8.2. Desincronización, `N=4`, `T=28`

| Indicador | Sin control `u=0` | Candidato de corrientes |
|---|---:|---:|
| Costo total | `92.0000` | **`30.238066`** (verificado con Tsit5) |
| Reducción relativa de costo | — | **`67.13 %`** |
| `r²(T)` | `≈1` | **`0.00063582`** |
| `E(T)` | `0` (sincronía mantenida) | `2.3768615` |
| `max |u_i|` | 0 | `0.8` |
| `max |Σ_i u_i|` | 0 | `≈1.11×10⁻¹⁶` |

El algoritmo utilizó varias semillas no nulas (ganadora `6`) para evitar aceptar automáticamente el extremal trivial `u=0`. Costo discreto RK4 `30.24694`; residuo máximo de ley PMP `≈1.05164`.

**Interpretación:** el candidato cumple las restricciones de amplitud y balance, reduce el costo y deja el orden de fase global muy bajo. Sin embargo, `r²(T)≈0` no descarta grupos de fases opuestas; tampoco se probó desincronización sostenida para cualquier horizonte, ni optimalidad global.

**Figuras:** [E_01: sincronización candidata](figures/E_01_sincronizacion_candidata.png), [E_02: residuos PMP de sincronización](figures/E_02_residuos_PMP_sync.png), [E_03: desincronización candidata](figures/E_03_desincronizacion_candidata.png), [E_04: residuos PMP de corrientes](figures/E_04_residuos_PMP_desync.png).

## 9. Etapa F — Hamilton–Jacobi–Bellman exploratoria

Las HJB completas dependen de `2N` estados. Se ensayaron dos problemas **distintos** en `N=2`; ninguno constituye una política global verificada para la red de `N=4` de la Etapa E.

### F1. Sincronización: reducción exacta especial en 2D

Con `c=(0,1.4)` y `N=2` sobre una subvariedad antipodal invariante, la dinámica de diferencias se reduce exactamente a dos variables `(d_v,d_w)`:

$$
\dot d_v=(1+\lambda)d_v-\frac{d_v^3}{3}-d_w,\qquad
\dot d_w=\frac{d_v-bd_w}{\tau},\qquad
E=d_v^2+d_w^2.
$$

Para estados antipodales no centrales se tiene **`r²=0`**, de modo que este caso no describe sincronización general de fase oscilatoria.

| Indicador | Aproximación denominada HJB | Referencia `λ=-1.5` |
|---|---:|---:|
| Llegada a `E≤0.03` | `0.968564` | `1.777938` |
| Costo reportado `J` | `3.410716` | `3.780949` |
| Celdas con función de valor finita | **59/2091 (2.82 %)** | No aplica |
| Valor aproximado en el estado inicial `d₀` | **`+∞`** | No aplica |
| Pasos fuera de cobertura de política | **14** | No aplica |

Aunque la simulación etiquetada «HJB» llegó antes y con menor costo, **no puede afirmarse que la mejora la produjo una realimentación HJB bien resuelta a lo largo de la trayectoria**: el valor es infinito en el estado inicial y hubo pasos fuera de cobertura. Estos números documentan el experimento, no una comparación concluyente de políticas óptimas.

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

- `videos/F_01_comparacion_control_variable_vs_constante.mp4`: mismo `N=4`, estado inicial y modelo; control variable factible versus `λ=-1.5`.
- `videos/F_02_comparacion_corrientes_vs_sin_control.mp4`: mismo `N=4`, estado inicial y modelo; candidato de corrientes versus `u=0`.

**No son comparaciones HJB vs PMP**: emplean las trayectorias de la Etapa E. Cuando la rama controlada llega antes a la tolerancia y se congela visualmente, ello no representa una continuación física real tras ese tiempo. Los MP4 no estaban incluidos en el paquete de revisión; sus nombres, existencia y tamaños constan en el inventario.

## 10. Resultados demostrados, calculados y pendientes

| Afirmación | Estado lógico | Precisión importante |
|---|---|---|
| Equilibrio individual único para parámetros admisibles | Analítico | Demostración en el cuaderno teórico |
| Condiciones locales de Hopf y fórmula de `κ` | Analítico | Dependen de las hipótesis para Hopf no degenerada |
| Existencia de alguna órbita periódica en régimen especificado | Analítico bajo hipótesis | No asegura unicidad ni estabilidad global |
| `λ≤−1−ε` con `b>0` implica decaimiento exponencial transversal | Analítico | Condición suficiente para igualdad de estados, no de actividad oscilatoria |
| Con corriente común acotada no se alcanza sincronía de estados exacta en tiempo finito desde fuera | Analítico bajo hipótesis | Justifica la tolerancia `δ>0` |
| Un control de corrientes diferenciadas rompe la igualdad inicial | Analítico y ensayos | No garantiza fases persistentemente desalineadas |
| Ciclo y multiplicadores de Floquet en `I=1.4` | Numérico | Shooting y validación por tests |
| Controles E más económicos que sus referencias | Numérico en casos concretos | Comparar solo problemas equivalentes; factibilidad verificada |
| Candidatos E satisfacen todas las condiciones PMP con residuo pequeño | **No demostrado** | Residuos cercanos a `1`; requiere refinamiento |
| Políticas F son soluciones globales válidas de HJB | **No demostrado** | Cobertura y convergencia insuficientes |
| Baja `r²` implica ausencia de clústeres | **Falso en general** | Clústeres en oposición pueden cancelar el orden global |
| La tolerancia de estados implica alineamiento de fases sin hipótesis adicionales | **No en general** | Ejemplo E: `E(T)≈0.03` con `r²(T)≈0.000185` |

## 11. Figuras y videos sugeridos para el informe/presentación

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

**Inventario completo del paquete de revisión:** 29 PNG (A: 3; B: 4; C: 6; D: 9; E: 4; F: 3) y 5 videos MP4 listados (C: 2; D: 1; F: 2). Los videos no fueron transferidos dentro del ZIP; deben abrirse desde `docs/videos/` para seleccionar fotogramas y evaluar velocidad/legibilidad.

## 12. Conclusiones aptas para el informe

1. El estudio individual identifica parámetros con dinámica oscilatoria estable y permite calcular críticamente puntos de Hopf, un ciclo límite y su estabilidad orbital.
2. El signo del acoplamiento determina un cambio cualitativo importante: valores negativos inducen atracción transversal, mientras los positivos estudiados sostienen separación entre neuronas. La cota global de Lyapunov es conservadora; se observó sincronización incluso cuando no es aplicable.
3. La desincronización por corrientes individuales es posible sin alterar el estímulo medio. Las estrategias manuales en grupos inducen separación de estados, pero no necesariamente un mejor costo ni desalineación final fuerte.
4. La optimización numérica produjo **candidatos admisibles y mejores que controles de referencia** en dos escenarios concretos. Son resultados cuantitativos relevantes, aunque no pruebas de optimalidad global o satisfacción precisa de PMP.
5. Las métricas de sincronización son complementarias: el error `E` mide proximidad absoluta de estados y `r` mide direcciones respecto de un centro. El experimento E de sincronización enseña que **una tolerancia pequeña de estados puede coexistir con `r` muy bajo**.
6. Las formulaciones HJB están matemáticamente establecidas y se implementaron prototipos reducidos; las simulaciones actuales muestran por qué la **cobertura del espacio de estados** y la maldición de dimensionalidad limitan su uso directo.

## 13. Limitaciones y trabajo futuro (orden sugerido)

1. **Interpretación físico-dinámica:** estudiar si las neuronas siguen oscilando al sincronizarse o convergen hacia el equilibrio; controlar distancias al centro geométrico y considerar una fase adaptada al ciclo límite.
2. **Sincronización conjunta estados–fases:** estudiar una restricción terminal adicional sobre `r` o un costo/criterio que evite resultados con `E` pequeño y fases geométricas opuestas, justificando primero su factibilidad.
3. **Pontryagin:** refinar resolución de adjuntas, transversalidad, tiempo final, multiplicadores y ley proyectada; comparar una transcripción más fina con shooting indirecto o métodos especializados.
4. **HJB:** experimentos de convergencia al refinar mallas, pasos y discretización de acciones; mejorar cobertura, interpolación y tratamiento de frontera antes de interpretar la realimentación obtenida.
5. **Clústeres y desincronización sostenida:** medir orden por subgrupos, distribución de diferencias de fase, permanencia temporal de `r²` bajo umbrales y efectos de controles no constantes.
6. **Robustez:** variar parámetros, tamaño y topología de la red, condiciones iniciales, tolerancias y posibles perturbaciones; repetir comparaciones con idénticos presupuestos de control.
7. **Comunicación:** inspeccionar los videos MP4 originales, seleccionar fotogramas, normalizar escalas y añadir anotaciones de `λ(t)`, `u(t)`, `E(t)` y `r(t)` a material de exposición.

## 14. Cómo reproducir / dónde están los datos

Ejecutar desde la raíz del proyecto, con el entorno Julia ya configurado:

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. tests/test_stage_C.jl
julia --project=. tests/test_stage_D.jl
julia --project=. tests/test_stage_E.jl
julia --project=. tests/test_stage_F.jl
```

Los experimentos originales se lanzaron con `scripts/run_stage_A.jl` a `scripts/run_stage_F.jl`; los videos comparativos desde `scripts/run_stage_F_video_only.jl`. Para consolidar resultados sin volver a optimizar:

```powershell
julia --project=. scripts/recopilar_resultados.jl
```

**Archivos fuente principales:**

- `data/processed/etapa_B_resultados.jld2` (órbita y Hopf).
- `data/processed/etapa_C_lambda_*.jld2` y `C_notebook_lambda_m1p5.jld2` (red).
- `data/processed/etapa_D_*.jld2` (controles y costos prefijados).
- `data/processed/etapa_E_sync.jld2`, `etapa_E_desync.jld2` (candidatos).
- `data/processed/etapa_F_hjb.jld2` (grillas HJB y simulaciones).
- `docs/figures/` (29 figuras según inventario); `docs/videos/` (5 MP4 según inventario).

**Fuente teórica primaria:** el enunciado original de MA4703 *Control de sincronización en redes de neuronas tipo FitzHugh–Nagumo* y el *Cuaderno teórico completo* (objetivos 1–5). El resumen ejecutivo del mismo proyecto complementa, pero no sustituye, las demostraciones del cuaderno.

**Revisión futura:** actualizar este Markdown únicamente con resultados calculados y contrastados; no transformar los pendientes en hechos concluidos. Se recomienda conservar el ZIP de revisión como una fotografía reproducible del estado del proyecto a fecha 2026-11-15.
