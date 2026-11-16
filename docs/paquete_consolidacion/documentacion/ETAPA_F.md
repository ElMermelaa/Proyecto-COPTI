# Etapa F — Hamilton–Jacobi–Bellman (MA4703)

## Alcance honesto

El cuaderno teórico formula dos ecuaciones:

**Sincronización:** `0 = min_{λ∈Λ} { L_s(x,λ)+∇W(x)·(F₀(x)+λG(x)) }`, con `W=0` en `Kδ={E≤δ}` (tiempo final libre).

**Desincronización:** `−V_t = min_{u∈U_I} { α_d r²(x)+γ_d||u||² +∇V(t,x)·(F_{λ₀}(x)+Bu) }`, con `V(T,x)=κr²(x)`.

Las dos HJB completas viven en dimensión `2N`. Aquí se hacen **dos experimentos exploratorios diferentes**, sin prometer resolver las EDP de alta dimensión ni probar optimalidad global.

## F1. Sincronización, reducción exacta a 2 dimensiones

Con los parámetros `a=0.7, b=0.5, τ=12.5, I=1.4`, el equilibrio único es `c=(0,1.4)`. Para `N=2` y corrientes comunes, la subvariedad

`z₁=(dᵥ,c𝓌+d𝓌), z₂=(−dᵥ,c𝓌−d𝓌)`

es **invariante** por simetría impar de `f(v)=v−v³/3`, ya que `cᵥ=0`.

La dinámica reducida exacta es

```
dᵥ' = (1+λ)dᵥ − dᵥ³/3 − d𝓌
d𝓌' = (dᵥ − b d𝓌)/τ
```

Aquí `E=dᵥ²+d𝓌²`, y para `d≠0` los vectores de fase son opuestos: **`r²=0`**. Por eso el costo fuera del objetivo es `ρ+αE+β+γ(λ−λref)²` y `W=0` si `E≤δ`. Para `δ>0` el proceso se detiene antes de llegar al centro singular.

**Limitación física decisiva:** esta reducción ilustra un mecanismo de llegada de estados hacia el equilibrio común. **No modela sincronización de oscilaciones con fases alineadas**; no puede sustituir los experimentos de `N=4` de la Etapa E. Una política HJB en esta subvariedad no es una realimentación óptima global para la red completa.

Se resuelve mediante iteración semi-lagrangiana, con `λ` muestreado finitamente y interpolación bilineal. Fuera de la grilla no se extrapola el valor; esas regiones aparecen como `Inf` y no se consideran resueltas.

## F2. Desincronización N=2, cuatro dimensiones completas

Se utiliza el campo exacto de la red `N=2`, conservando las cuatro variables `(v₁,w₁,v₂,w₂)`. Dado `u₁+u₂=0`, se escribe `u=(q,−q)` y se busca `q∈[−Umax,Umax]`.

El problema discreto implementa Bellman hacia atrás, sobre una grilla 4D y horizonte corto (`T=1.2` por defecto):

```
V_k(x) ≈ min_q { h [α_d r²(x)+2γ_d q²] + interp(V_{k+1}, Φ_h(x,q)) }
V_M(x) = κr²(x).
```

`Φ_h` aproxima el flujo del sistema mediante RK2. Se excluyen fases indefinidas y transiciones fuera de la malla. No se extrapola fuera del dominio. Las corrientes son **siempre** balanceadas. El problema completo de Etapa E usa `N=4`, `T=28`: no se debe afirmar que el control de F2 lo resuelva ni comparar directamente sus costos.

Desde una condición exactamente sincronizada pueden existir simetrías, múltiples minimizadores y funciones de valor no suaves. Si la realimentación numérica devuelve `u=0`, no implica que sea imposible desincronizar: la discretización también puede no detectar un rompimiento de simetría ventajoso.

## Comparación audiovisual

`run_stage_F.jl --video` usa **la caché de sincronización de la Etapa E** (`data/processed/etapa_E_sync.jld2`) y genera una comparación de dos controles sobre la **misma red N=4 y el mismo estado inicial**:

- Control variable factible calculado en E.
- Control constante `λ=-1.5`.

Además, si existe `data/processed/etapa_E_desync.jld2`, genera una segunda comparación de **corrientes controladas frente a `u=0`** en la misma red `N=4`, con el mismo acoplamiento y estado inicial.

El control variable se detiene al llegar a su `T` y su estado se **congela solo para la visualización**. No se simula una extensión ficticia del control. En esta animación no se incluye una política HJB de F1 porque **no representa el mismo experimento**.

## Ejecución y validación

Desde la raíz del repositorio:

```
julia --project=. tests/test_stage_F.jl
julia --project=. scripts/run_stage_F.jl
julia --project=. scripts/run_stage_F_video_only.jl
```

Si falta `etapa_E_sync.jld2`, primero corre `scripts/run_stage_E.jl`. `run_stage_F_video_only.jl` carga los resultados E y NO recalcula las grillas HJB. El argumento `--video` en `run_stage_F.jl` es una alternativa que sí vuelve a calcular las grillas.

El test cubre reducción matemática, construcción de valores, restricciones de corrientes y salidas finitas. No reemplaza un estudio de convergencia de malla, estabilidad/consistencia del esquema, validación de soluciones de viscosidad ni verificación de optimalidad. Para mejorar precisión, comparar distintas grillas `vgrid/wgrid`, `vaxis/waxis`, acciones y pasos temporales.
