# Etapa E — Control óptimo, Pontryagin y controles candidatos (MA4703)

## 1. Alcance y decisiones

Se conservan **sin cambios** el modelo del enunciado, la convención de acoplamiento y los dos problemas del cuaderno teórico. Esta entrega amplía la Etapa D sin sustituir sus archivos analíticos ni sus resultados guardados.

Modelo en una red completa idéntica:

\[
\dot v_i=v_i-v_i^3/3-w_i+I_0+u_i+\lambda(v_i-\bar v),\qquad
\dot w_i=(v_i+a-bw_i)/\tau.
\]

Se usan el estado intercalado `x=[v1,w1,v2,w2,...]`, centro de fase fijo `c`, error de estados `E=(1/N)Σ||zi-z̄||²` y `r²=||(1/N)Σ(zi-c)/||zi-c||||²` sólo donde todas las fases están definidas. En particular, la penalización de sincronización es **β(1−r²)²**.

### A. Sincronización — tiempo libre

\[
\min_{T>0,\ \lambda(t)\in[\lambda_{min},\lambda_{max}]}
\int_0^T [\rho+\alpha E(x)+\beta(1-r(x)^2)^2+
\gamma(\lambda-\lambda_{ref})^2]dt
\]

sujeto a `x(0)=x0` y `E(x(T))≤δ`, con `0<δ<E(x0)`. La restricción debería ser activa en un óptimo verdadero de tiempo libre. El control normal necesario es

\[
\lambda^*(t)=\Pi_{[\lambda_{min},\lambda_{max}]}
\left(\lambda_{ref}-\frac{1}{2\gamma}\sum_i p_{v_i}(v_i-\bar v)\right).
\]

Transversalidad normal: `p(T)=μ∇E(x(T))`, `μ≥0`, `E(T)=δ`, `H(T)=0`. **Estas son condiciones necesarias** bajo hipótesis de existencia, normalidad y regularidad; no constituyen una solución numérica por sí mismas.

### B. Desincronización — tiempo fijo

\[
\min_{u(\cdot)}\int_0^T[\alpha_d r^2(x(t))+
\gamma_d\sum_i u_i(t)^2]dt+\kappa r^2(x(T)),
\quad |u_i|\le U_{max},\ \sum_i u_i=0.
\]

`λ=λ0` es fijo. La ley PMP necesaria minimiza un problema cuadrático estrictamente convexo **por instante**:

\[
u_i^* = \Pi_{[-U_{max},U_{max}]}\big(-(p_{v_i}+\nu)/(2\gamma_d)\big),
\quad \sum_i u_i^*=0.
\]

Se calcula `ν` mediante bisección, **no** se recortan las corrientes independientemente después de centrarlas. La condición terminal es `p(T)=κ∇r²(x(T))`. Desde sincronía exacta el control `u=0` es un extremal estacionario, por lo que se exploran *varias semillas no nulas* y se compara siempre con `u=0`.

## 2. Método numérico implementado

1. **Transcripción directa reducida:** el control es constante en `M` subintervalos. Sincronización también optimiza `T`; desincronización mantiene `T` fijo. El campo y el costo acumulado se integran juntos mediante RK4 con `substeps` por intervalo.
2. **Derivadas automáticas:** `ForwardDiff.gradient` de la función objetivo discretizada respecto de todos los tramos (y `T` en sincronización).
3. **Restricciones:** proyección de λ y T sobre intervalos; desincronización proyecta cada vector de corrientes sobre la caja intersectada con el hiperplano `Σu=0`.
4. **Sincronización:** Lagrangiano aumentado de desigualdad `E(T)-δ≤0`, con multiplicador aproximado `μ≥0`. Después de cada minimización interna se actualiza el multiplicador y la penalización. **Puede terminar sin factibilidad**: verificar `feasible`, el residuo terminal y el error de integración independiente.
5. **Desincronización:** descenso proyectado multistart con seis semillas deterministas no nulas (menos para N=2), más `u=0`. Devuelve el mejor costo discreto obtenido, incluso si sigue siendo el control nulo. Un candidato no trivial **no** implica `r(T)≈0`.
6. **Verificación externa:** los candidatos se integran otra vez con `OrdinaryDiffEq.Tsit5()` usando `tstops` en los saltos. Los resultados se comparan con RK4. Se integran **hacia atrás las ecuaciones adjuntas continuas** sobre la trayectoria independiente; se informa el residuo de la ley minimizadora PMP y, para sincronización, `H(T)`.
7. **Shooting indirecto opcional:** `shoot_sync_pmp` resuelve un sistema de `2N+2` ecuaciones con incógnitas `p(0),T,μ`. Se recomienda probar N=2 y buenas semillas. El solver puede no converger: se registra su estado; nunca se declara automáticamente óptimo un cero numérico.

No se utilizaron `OptimalControl.jl`/`JuMP.jl` todavía porque los primeros experimentos funcionan con las dependencias existentes. Una comparación con transcripción directa de orden superior y optimizadores especializados es una ampliación natural. Para problemas grandes, el descenso proyectado implementado aquí es **un prototipo docente**, no una biblioteca de optimización de producción.

## 3. Estructura de archivos

```text
src/FitzHughNagumoControl.jl        (único archivo existente a reemplazar)
src/control/optimal_pmp.jl           (Hamiltonianos, derivadas, controles PMP)
src/control/optimal_direct.jl        (transcripción RK4 y optimizadores)
src/analysis/optimal_diagnostics.jl  (Tsit5 independiente, adjuntas, shooting)
src/visualization/optimal_plots.jl  (figuras de J, E, r², controles y PMP)
src/utils/optimal_storage.jl        (guardar en JLD2 sin closures)
tests/test_stage_E.jl               (pruebas matemáticas y numéricas)
scripts/run_stage_E.jl              (experimento reproducible)
notebooks/05_pontryagin.ipynb       (análisis interactivo mediante celdas)
```

## 4. Cómo ejecutar

Desde la raíz del repositorio, con `Project.toml` de las etapas anteriores:

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. tests/test_stage_C.jl
julia --project=. tests/test_stage_D.jl
julia --project=. tests/test_stage_E.jl
julia --project=. scripts/run_stage_E.jl
```

La Etapa E requiere `OrdinaryDiffEq`, `NonlinearSolve`, `ForwardDiff`, `CairoMakie`, `JLD2`, que ya se instalaron en etapas previas. No reescribir `Project.toml`/`Manifest.toml`.

Resultados principales: `docs/figures/E_01...E_04*.png`, `data/processed/etapa_E_sync.jld2` y `etapa_E_desync.jld2`. Ejecutar `--shooting` sólo si se desea intentar el método indirecto adicional.

## 5. Cómo interpretar resultados

- `result.J`: costo RK4 del control encontrado.
- `result.feasible`: indicador numérico de que `E(T)-δ` está por debajo de la tolerancia elegida. Es necesario **verificar con Tsit5** el valor real `E(T)` y si hubo llegada antes de T.
- `result.mu`: multiplicador aproximado de la desigualdad terminal (no necesariamente convergido).
- `diag.max_projected_control_residual`: diferencia máxima entre control aplicado y ley PMP calculada desde adjuntos reconstruidos (muestreada).
- `diag.terminal_H`: error de la condición `H(T)=0` de sincronización con tiempo libre.
- `result.baseline_J` de desincronización: costo del control `u=0`, con exactamente los mismos pesos y horizonte.
- `r²≈0` indica cancelación vectorial de fases, **no** ausencia de clústeres opuestos; `E≈0` indica estados próximos, **no** permanencia oscilatoria.

Una buena aproximación numérica a una condición necesaria no equivale a prueba de optimalidad global. La HJB y las verificaciones analíticas se abordan en la Etapa F.

## 6. Limitaciones y tareas futuras

- El gradiente de fase requiere `||zi-c||>center_tol`; una trayectoria que toque el centro no es admisible para el costo geométrico tal como fue definido.
- Los controles son constantes a trozos; el descenso proyectado puede encontrar mínimos locales y estancarse con restricciones activas.
- El Lagrangiano aumentado no garantiza satisfacción de la condición terminal dentro de un número prefijado de iteraciones; se reporta si falla.
- El extremo trivial `u=0` puede ser seleccionado si todos los controles ensayados cuestan más.
- La dinámica de la red sincronizada puede converger a un equilibrio; no se está imponiendo preservación de ciclos límite.
- Aumentar `M` y `substeps`, reducir tolerancias de integración, comparar valores del costo y repetir con diversas semillas es necesario para realizar un estudio de convergencia serio.
