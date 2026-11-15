# Etapa B — Bifurcaciones de Hopf y ciclos límite en Julia

## 1. Modelo y dependencia de la Etapa A

Mantenemos exactamente las ecuaciones del curso:

\[
\dot v=v-v^3/3-w+I,\qquad \dot w=(v+a-bw)/\tau.
\]

Esta entrega **extiende** la Etapa A: no modifica las ecuaciones, parámetros ni
el código de simulación validado con 61 pruebas. Sí actualiza el archivo central
`src/FitzHughNagumoControl.jl` para incluir/exportar módulos B. **No sobreescribir**
`Project.toml` ni `Manifest.toml` existentes.

## 2. Bifurcaciones y ecuaciones verificables

Para `0<b<1` y `tau>b`, con `s=sqrt(1-b/tau)`:

\[
v_H^\pm=\pm s,\quad
I_H^\pm=(v_H^\pm)^3/3-v_H^\pm+(v_H^\pm+a)/b,\quad
\omega_H=\sqrt{\tau-b^2}/\tau.
\]

El cruce transversal cumple

\[
\frac{d\sigma}{dI}\bigg|_{I_H}=
-\frac{b v_H}{1-b+bv_H^2}\ne0,
\quad\sigma=\frac12\operatorname{tr}J^*.
\]

Según el cuaderno teórico (capítulo 2):

\[
\kappa=\frac{(2b-1)\tau-b^2}{8(\tau-b^2)}.
\]

- `kappa < 0`: Hopf **supercrítica**, ciclo localmente estable del lado de equilibrio inestable.
- `kappa > 0`: Hopf **subcrítica**, ciclo pequeño inestable del lado de equilibrio estable.
- `kappa = 0`: coeficiente degenerado; no concluir usando solo la forma normal cúbica.

Los valores de `hopf_currents` / `hopf_analysis` son candidatos analíticos de
Hopf, no resultados numéricos de continuación.

Para los parámetros del ejemplo de la Etapa A `a=0.7, b=0.5, tau=12.5`:

- `I_minus ≈ 0.1066694158`, `I_plus ≈ 2.6933305842`.
- `omega_H = 0.28` y `T_lineal = 2pi/0.28 ≈ 22.44`.
- `kappa ≈ -0.00255102`: Hopf supercrítica.

El período lineal en Hopf **no** es necesariamente el período del ciclo de gran
amplitud para `I=1.4`.

## 3. ¿Cómo detectamos ciclos? (dos niveles)

### Primer nivel: cruces de Poincaré y amplitudes

- Integrar hasta `transient+window`; descartar el transitorio.
- Usar la sección `v=v*(I)` y registrar cruces ascendentes mediante interpolación densa.
- Calcular los últimos tres períodos de retorno y las amplitudes `max(v)-min(v)`.
- Verificar que ambos sean suficientemente estables.

Se devuelve `status=:periodic_candidate` solamente si hay suficientes cruces,
la amplitud no es despreciable y los períodos y amplitudes se estabilizan.
Una `:insufficient_crossings` o `:not_converged` **no demuestra** inexistencia.

### Segundo nivel: shooting + sistema variacional

Se resuelve el sistema de tres ecuaciones

\[
\Phi_T(z_0)-z_0=0,\qquad (z_0)_1-v^*=0,
\]

para las tres incógnitas `v(0),w(0),T`. La sección fija la fase y evita
indeterminación por traslación temporal.

El Jacobiano de Newton utiliza la matriz fundamental

\[
\dot M=DF(z(t))M,\qquad M(0)=I_2,
\]

con derivada final de flujo respecto al período igual a `F(z(T))`.
El resultado numérico reporta residuo de cierre, amplitudes y multiplicadores
de Floquet. El multiplicador tangencial debe ser aproximadamente uno; el otro
permite clasificar la estabilidad orbital. Un resultado `:converged` es un
ciclo numérico corregido, **no** un teorema de unicidad/globalidad.

## 4. Continuación de rama y límites

`continue_periodic_branch(p, orbit, currents)` predice una órbita próxima en
cada nueva corriente y la corrige por shooting, usando el punto anterior.
Sirve para ramas accesibles en pasos suficientemente pequeños; **no** es
pseudo-arclength y puede fallar en pliegues, zonas rígidas, o cerca de órbitas
muy pequeñas. Un fallo se reporta en `status` y detiene la rama.

`hopf_seed(p, :minus/:plus; offset=...)` devuelve una aproximación local de
forma normal que puede servir como inicialización de *shooting*, incluso para
ciclos inestables. No garantiza que Newton encuentre el ciclo. Para realizar
continuación robusta de ramas inestables, pliegues y diagramas completos se
recomienda posteriormente `BifurcationKit.jl` (shooting o collocation con
pseudo-arclength), que aquí no se instala para evitar una dependencia grande.

## 5. Estructura de archivos de esta entrega

```
src/FitzHughNagumoControl.jl                 # actualizado
src/analysis/bifurcations.jl                # NUEVO
src/visualization/bifurcation_plots.jl      # NUEVO
scripts/run_stage_B.jl                      # NUEVO
tests/test_stage_B.jl                       # NUEVO
notebooks/02_bifurcations.ipynb             # NUEVO
docs/ETAPA_B.md                             # NUEVO
```

Además, `data/processed/etapa_B_resultados.jld2` y cuatro figuras
`docs/figures/B_*.png` se generan al ejecutar el script. Son **resultados**,
no archivos de código entregados, y no se han generado aquí con Julia.

## 6. Instalación y ejecución en Windows (raíz del repositorio)

Añade únicamente una dependencia nueva (para guardar resultados). En **CMD o PowerShell**,
abre primero Julia con:

```powershell
julia --project=.
```

Dentro de `julia>` ejecuta:

```julia
using Pkg
Pkg.add("JLD2")
exit()
```

**No** escribas `julia --project=.` dentro de `julia>`: ese comando es de la
terminal de Windows.

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. scripts/run_stage_B.jl
```

Después abre `notebooks/02_bifurcations.ipynb`, selecciona el kernel Julia
y **reinícialo** si anteriormente habías cargado módulos A.

## 7. Qué resultados deberían observarse

1. Cruces de traza en las dos corrientes analíticas de Hopf.
2. Para `b=0.5, tau=12.5`, κ negativa (supercrítica).
3. Para `I=1.4`, tras el transitorio, una oscilación estable con período
   medible y cierre por shooting.
4. Un multiplicador Floquet cercano a uno y otro de módulo menor que uno.
5. Un barrido que muestre candidatos de ciclo solo donde la detección fue
   satisfactoria. Los valores no detectados se conservan como `NaN`.
6. Una continuación corta de la órbita hacia corrientes cercanas; sus registros
   pueden mostrar fallos en lugar de ocultarlos.

**Precisión científica:** la existencia de al menos una órbita periódica bajo
las hipótesis probadas se apoya en Poincaré–Bendixson (cuaderno); su unicidad y
estabilidad global no están demostradas. Los cálculos numéricos no alteran ese
alcance.

## 8. Fuentes

- Cuaderno teórico completo del proyecto MA4703: capítulos 1–2, enunciado original.
- OrdinaryDiffEq: https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/usage/
- Referencia de métodos de órbitas periódicas de BifurcationKit:
  https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/periodicOrbit/
- JLD2: https://juliaio.github.io/JLD2.jl/stable/basic_usage/
