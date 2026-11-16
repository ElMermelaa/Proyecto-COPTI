# Etapa D — Controles suficientes, cotas de Lyapunov y evaluación de costos

Proyecto MA4703 · Universidad de Chile · Julia

## 1. Alcance científico

La Etapa D implementa el **objetivo 2**, sin afirmar optimalidad. Reutiliza:

- Modelo de la Etapa A: `v̇=v-v³/3-w+I`, `τẇ=v+a-bw`.
- Bifurcaciones y ciclos refinados de la Etapa B.
- Redes, métricas y visualización de la Etapa C.

En todos los módulos, el signo **se conserva**:

`v̇_i=f(v_i)-w_i+I0+u_i+λ(v_i-v̄)`.

Un acoplamiento `λ<0` es atractivo.

## 2. Dos leyes suficientes de sincronización

Definimos `e_i=v_i-v̄`, `h_i=w_i-w̄` y

`V=1/2 Σe_i² + τ/2 Σh_i²`, `E=(1/N)Σ(e_i²+h_i²)`.

La demostración del cuaderno establece, con corriente basal común y conectividad
exactamente todos-con-todos:

`V̇ ≤ (1+λ) Σe_i² - b Σh_i²`.

### 2.1 Acoplamiento constante

Si `b>0`, `τ>0` y `λ<-1`, elegimos `ε=-1-λ>0`:

`c=min(ε,b/τ)>0`, `V(t)≤V(0)exp(-2ct)`.

### 2.2 Corrientes diferenciadas en realimentación

`u_i(x)=-k(v_i-v̄)` con `k>1+λ0`, y sin saturación:

`V̇≤(1+λ0-k)Σe_i²-bΣh_i²`, `c=min(k-1-λ0,b/τ)`.

**Equivalencia de campos**: si `A_ij=1`, la dinámica con `λ0` y
`u_i=-k e_i` coincide con la dinámica sin esas corrientes y con
**`λ_efectivo=λ0-k`**. Esto NO implica equivalencia de los costos de control,
ni se extiende sin más a una matriz de conectividad arbitraria.

Las garantías de esta implementación requieren `b>0`, corriente basal
constante y común, y `A_ij=1` fuera de la diagonal. El acoplamiento puede
ser constante o de la familia `AttractiveModulatedCoupling(ceiling,depth,omega)`:
`λ(t)=ceiling-depth*sin(omega*t)^2≤ceiling`; la cota es **analítica**.
No se extiende automáticamente a cualquier función `lambda(t)`.
Las métricas pueden estudiarse sin tales condiciones, **pero no se emite
certificado** sin ellas. El código no infiere cotas globales a partir de un
muestreo temporal de un control arbitrario.

### Tiempo suficiente

Para `0<δ<E(0)` definimos `Cτ=max(1,τ)/min(1,τ)` y

`T_suf=max{0, ln(Cτ E(0)/δ)/(2c)}`.

Este es un límite conservador: **no** predice el tiempo óptimo y **no**
certifica que las neuronas continúen sobre un ciclo límite.

## 3. Ruptura de sincronía

Para `λ(t)=λ0` usamos un patrón balanceado `Σu_i(t)=0`, construido con
`balanced_current_profile` y `balanced_current_law`.

En sincronía exacta, `v_i=v_j`, `w_i=w_j`, se tiene

`d(v_i-v_j)/dt = u_i-u_j`.

Por tanto una corriente diferenciada rompe la igualdad de estados, aunque
**no implica automáticamente** `r(t)≈0` ni fases desincronizadas de manera
sostenida. El experimento empieza exactamente en la variedad sincronizada
(y, cuando está disponible, con un punto de un ciclo refinado de la Etapa B).

Para las corrientes balanceadas elegidas, `|u_i|≤amplitude`. No se aplica
feedback saturado aquí: saturar entradas individualmente después de centrar
puede destruir `Σu_i=0`.

## 4. Funciones objetivo: solo evaluación, no minimización

La Etapa D añade código reutilizable para las formulaciones adoptadas del
cuaderno (objetivo 3), pensando en la Etapa E.

### Sincronización (tiempo libre)

`J_s(T,λ)=∫₀ᵀ [ρ + αE + β(1-r²)² + γ(λ-λref)²]dt`,

sujeto a `E(x(T))≤δ`.

`evaluate_sync_cost` calcula:

- `horizon_J`: integral hasta el final **prefijado de la simulación**; se usa
  solo como diagnóstico, no como costo admisible del problema de tiempo libre.
- `arrival_J`: integral hasta el **primer cruce detectado** de `E≤δ`.
  Vale `Inf` si no se detectó llegada en la simulación.
- `arrival_time` y `status` para no confundir `no alcanzado` con `costo bajo`.

El cruce se refina mediante bisección usando `sol(t)`; pueden perderse cruces
muy breves entre dos tiempos guardados. Usar `saveat` pequeño. Si alguna
neurona coincide con el centro de fase elegido, se lanza `DomainError`: el
costo del cuaderno **no está definido** en ese caso.

### Desincronización (tiempo fijo)

`J_d(u)=∫₀ᵀ [α_d r²+γ_d Σu_i²]dt + κ r²(T)`.

El evaluador comprueba en los tiempos muestreados `Σu_i≈0` y `|u_i|≤Umax`.
Para leyes arbitrarias dadas por funciones, este chequeo **no garantiza**
el cumplimiento entre muestras. La ley balanceada suministrada satisface
las cotas por construcción.

Los pesos `ρ,α,β,γ,λref,α_d,γ_d,κ` de las demostraciones son
**ejemplos numéricos configurables**, no datos fijados por el enunciado.

**No debemos comparar directamente** `J_s` con el esfuerzo de un feedback
por corrientes: `J_s` penaliza `λ`, no los `u_i` de la realimentación.
`total_current_effort` se reporta de forma separada. Tampoco comparar
`J_s` (sincronizar) con `J_d` (desincronizar) como si optimizaran lo mismo.

## 5. Experimentos generados

| Archivo PNG | Descripción |
|---|---|
| `D_01_certificado_acoplamiento.png` | V y E frente a sus cotas (λ=-1.5) |
| `D_02_certificado_feedback.png` | Cotas para `λ0=0, k=1.5` |
| `D_03_comparacion_sincronizacion.png` | Casos λ=-1.5,-0.3,-0.15,+0.45, λ(t) modulada y feedback |
| `D_04_corrientes_feedback.png` | Corrientes reales del feedback |
| `D_05_costos_J_sincronizacion.png` | J acumulado y cruces de tolerancia |
| `D_06_ruptura_sincronia.png` | Corrientes diferenciadas desde sincronía exacta |
| `D_07_corrientes_desincronizacion.png` | Señales u_i(t), λ0 |
| `D_08_costos_J_desincronizacion.png` | Integrales acumuladas de J_d |
| `D_09_grafo_con_corrientes.png` | Vista final de la red estimulada |

Video opcional: `docs/videos/D_01_corrientes_diferenciadas.mp4`.

Datos numéricos: `data/processed/etapa_D_*.jld2` (parámetros, tiempos,
estados, λ, corrientes, E, V, r, sumarios de J). Los cierres y el solver
**no se serializan**. Las trazas se pueden volver a graficar sin reintegrar.

## 6. Instalación y ejecución

1. Copiar el **contenido** del ZIP a la raíz de tu repositorio y combinar
   `src`, `docs`, `tests`, `notebooks` y `scripts`.
2. Reemplazar **solo** `src/FitzHughNagumoControl.jl` cuando Windows pregunte:
   añade los `include`/`export` de la Etapa D. Los módulos originales A-C
   no están incluidos ni se modifican.
3. Conservar `Project.toml`, `Manifest.toml` y resultados previos.
4. Ejecutar en PowerShell o CMD desde la raíz:

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. tests/test_stage_C.jl
julia --project=. tests/test_stage_D.jl
julia --project=. scripts/run_stage_D.jl
julia --project=. scripts/run_stage_D.jl --video
```

El video es opcional y tarda más por el renderizado. No se necesitan
bibliotecas nuevas más allá de las ya instaladas hasta la Etapa C.

Para el notebook, reiniciar kernel antes de ejecutar desde el comienzo,
seleccionar el kernel Julia asociado al `Project.toml` del repositorio y
abrir `notebooks/04_sufficient_controls.ipynb`. El notebook usa nombres de
módulo explícitos (`FHN.`) para evitar conflictos de importación tras
varios `include` en el mismo kernel.

## 7. Próximo paso

La Etapa E reutilizará los costos aquí implementados y añadirá controles
admisibles optimizados, las ecuaciones adjuntas y la comprobación de
condiciones de Pontryagin. Las leyes suficientes de la Etapa D serán
**baselines**, no supuestos óptimos.
