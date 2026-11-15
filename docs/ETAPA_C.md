# Etapa C — Redes FitzHugh–Nagumo (MA4703)

## 1. Modelo y signos

Para neuronas idénticas, con `x = [v1,w1,...,vN,wN]`, implementamos:

\[
\dot v_i=v_i-\frac{v_i^3}{3}-w_i+I_0+u_i(t)
+\frac{\lambda(t)}{N}\sum_{j=1}^{N}A_{ij}(v_i-v_j),
\qquad \dot w_i=\frac{v_i+a-bw_i}{\tau}.
\]

Con `A_ij=1`, el sumatorio dividido por `N` es `v_i-v̄`. Este es **el signo del enunciado**: `λ < 0` es atractivo. Por defecto `u_i=0`, `λ=-0.5` y la matriz es `ones(N,N)` (incluye bucles propios nulos en la ecuación). Un anillo es solo una **extensión experimental**.

**Hipótesis**: `a≥0`, `0≤b<1`, `τ>0`; ciertas garantías de sincronización requerirán `b>0` y corriente común.

## 2. Invariancia de la variedad sincronizada

La variedad `S={x: z1=⋯=zN}` es invariante si todas las neuronas reciben la misma corriente. Sobre ella se anula el acoplamiento para **cualquier** valor de `λ(t)`: todas siguen la EDO individual. La prueba numérica comprueba derivadas iguales y que las trayectorias de N neuronas coinciden con una solución individual dentro de la precisión numérica.

Con corrientes individuales distintas, incluso sobre `S` se tiene `d(v_i-v_j)/dt = u_i-u_j`, lo que puede romper instantáneamente igualdad de estados. No se confunde esto con un resultado de desincronización sostenida de fases.

## 3. Métricas y centro de fase

\[
E(x)=\frac1N\sum_i[(v_i-\bar v)^2+(w_i-\bar w)^2],
\qquad
V(x)=\frac12\sum_i(v_i-\bar v)^2+\frac{\tau}{2}\sum_i(w_i-\bar w)^2.
\]

Para un centro fijo `c`, definimos `q_i=(z_i-c)/||z_i-c||`, `r²=||(1/N)Σq_i||²` y `r=√r²`. El valor predeterminado del centro es el equilibrio individual `equilibrium(p)`. Este es un **centro geométrico elegido**, no una fase isócrona de un ciclo límite; para interpretaciones precisas hay que comprobar el régimen oscilatorio.

Si `||z_i-c|| <= center_tol` para cualquier neurona, la fase se marca como **indefinida (`NaN`)**, sin sustituirla silenciosamente por 0 o 1. Además, `r≈0` puede indicar dos clústeres en oposición.

El cuaderno prueba para corriente común y campo medio:

\[
\dot V\le (1+\lambda)\sum_i (v_i-\bar v)^2-b\sum_i(w_i-\bar w)^2.
\]

Si `b>0` y `λ(t)≤-1-ε`, entonces `V(t)≤V(0)e^{-2ct}` con `c=min(ε,b/τ)`. Las pruebas del código comprueban esa cota para un caso con `λ=-1.35`. No es un umbral necesario; valores menos negativos pueden sincronizar.

## 4. Estabilidad longitudinal, transversal y Floquet

En equilibrio homogéneo y campo medio se obtienen

\[
J_{\parallel}=\begin{pmatrix}1-v_*^2&-1\\1/\tau&-b/\tau\end{pmatrix},\qquad
J_{\perp}=\begin{pmatrix}1-v_*^2+\lambda&-1\\1/\tau&-b/\tau\end{pmatrix}.
\]

`mean_field_linear_blocks` informa estabilidad transversal local del equilibrio y estabilidad del sistema completo (que exige estabilidad longitudinal **y** transversal). Un equilibrio longitudinalmente inestable puede coexistir con sincronización de oscilaciones.

Si existe un ciclo refinado `z*(t)` de período `P` de la Etapa B, `transverse_floquet` integra la ecuación variacional con `1-v*(t)^2+λ` sobre un período. Los multiplicadores transversales se obtienen de la monodromía. Si ambos están estrictamente dentro del disco unidad, el **ciclo sincronizado** es transversalmente atractivo de manera local.

## 5. Visualizaciones

- `plot_network_snapshot`: grafo circular; color según potencial. Aristas semitransparentes para evitar saturación visual.
- `plot_network_timeseries`: series temporales de potencial y recuperación.
- `plot_network_phase`: trayectorias individuales en el plano `(v,w)`.
- `plot_network_metrics`: evoluciones `E,V,r`, con `NaN` cuando falta fase.
- `plot_coupling_comparison`: varias intensidades `λ` con el **mismo estado inicial**.
- `animate_network_dashboard`: animación MP4/GIF, grafo de nodos parpadeantes + plano de fase + `E(t)` + `r(t)`; las posiciones circulares y tamaños son **codificaciones visuales**, no magnitudes físicas.

En animaciones, color y tamaño cambian con `v_i(t)`: al alcanzar estados similares parpadean conjuntamente. Las aristas no representan transmisión neuronal con retardos ni actividad real medida; solo la conectividad del modelo.

## 6. Reproducibilidad y caché

El script intenta obtener fases iniciales equiespaciadas de `data/processed/etapa_B_resultados.jld2` (generado por la Etapa B). Solo usa ese archivo si `(a,b,τ,I)` coinciden. Si no está disponible, usa `distributed_initial_state`, una distribución determinista que **no** se presenta como ciclo límite.

Las simulaciones nuevas guardan tiempos, estados, metadatos, indicadores `E,V,r²,r`, centro geométrico, parámetros del solver, conectividad y descripción de los controles con `JLD2`. No serializamos `ODESolution` ni funciones Julia. Las figuras/vídeos pueden generarse desde la solución calculada, y los arrays guardados permiten reconstruirlos sin integrar otra vez (reconstrucción de gráficos con los datos, no del objeto solución original).

## 7. Cómo ejecutar

En terminal PowerShell/CMD **fuera** de `julia>` y dentro de la raíz del repositorio:

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. tests/test_stage_C.jl
julia --project=. scripts/run_stage_C.jl
julia --project=. scripts/run_stage_C.jl --video
```

El último comando es opcional y puede demorar más (codificación de 150 fotogramas). Para GIF desde Jupyter cambia `CREAR_GIF = true` en `notebooks/03_networks.ipynb`.

## 8. Archivos y compatibilidad

Solo se **reemplaza** `src/FitzHughNagumoControl.jl`, para registrar las nuevas funciones. Los demás archivos son **nuevos**. No se tocan `Project.toml`, `Manifest.toml`, pruebas A/B ni módulos previos.

Dependencias nuevas: **ninguna** si ya completaste las etapas A y B (`OrdinaryDiffEq`, `CairoMakie`, `ForwardDiff`, `NonlinearSolve`, `IJulia`, `JLD2`).

**Limitaciones:** no se afirma sincronización sostenida, formación de clústeres ni desincronización únicamente mirando la animación. Para topologías distintas de campo medio, las desigualdades demostradas en el PDF no se trasladan sin una nueva prueba. La Etapa D estudiará de manera sistemática controles suficientes y comparación de métricas.
