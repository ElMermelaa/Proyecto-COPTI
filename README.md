# MA4703 — Redes FitzHugh–Nagumo y control óptimo

## Etapa A: dinámica individual

Implementación Julia del sistema individual del enunciado MA4703:

\[
\dot v = v-\frac{v^3}{3}-w+I,\qquad
\dot w = \frac{v+a-bw}{\tau}.
\]

Los módulos, el notebook y las pruebas **no alteran** el modelo ni los signos del enunciado. Los valores predeterminados `a=0.7`, `b=0.5`, `tau=12.5`, `I=1.4` son un ejemplo editable y **no** provienen de los PDF.

### Integración con tu repositorio existente

Este ZIP contiene rutas relativas a la carpeta raíz del repositorio. **Copia su contenido dentro del repositorio, fusionando carpetas**, pero **no reemplaces** tus actuales `Project.toml`, `Manifest.toml`, `.gitignore`, otros notebooks ni ajustes de VS Code. La entrega no incluye esos archivos.

Si ya existe un `README.md` personalizado, conserva el tuyo y copia a él estas instrucciones o consulta `docs/ETAPA_A.md`.

### Dependencias

Se aprovecha el entorno Julia que ya instalaste. Para Etapa A se usan:

- `OrdinaryDiffEq` para integrar las EDO.
- `CairoMakie` para los gráficos.
- `ForwardDiff` para verificar el jacobiano en los tests.
- `LinearAlgebra` y `Test`, bibliotecas estándar de Julia.
- `IJulia` para ejecutar el notebook.

`NonlinearSolve` puede estar instalado, pero en esta etapa la **bisección monotónica** resuelve el equilibrio sin depender de un algoritmo no lineal externo. Lo emplearemos más adelante cuando sea útil.

### Ejecución rápida

Desde la terminal de VS Code ubicada en la **raíz de tu repositorio**:

```powershell
julia --project=. tests/runtests.jl
julia --project=. scripts/run_stage_A.jl
```

El segundo comando genera:

- `docs/figures/A_01_series_temporales.png`
- `docs/figures/A_02_retrato_fase.png`
- `docs/figures/A_03_estabilidad_vs_corriente.png`

Después abre `notebooks/01_individual.ipynb` en VS Code, elige el kernel de Julia y ejecuta las celdas en orden. El notebook detecta automáticamente la raíz del repositorio buscando `src/FitzHughNagumoControl.jl`.

### Organización

| Ruta | Función |
|---|---|
| `src/FitzHughNagumoControl.jl` | módulo principal y exportaciones |
| `src/models/parameters.jl` | parámetros y validación |
| `src/models/individual.jl` | campo, nulclinas y solucionador |
| `src/analysis/equilibria.jl` | equilibrio único y residuo |
| `src/analysis/stability.jl` | jacobiano, autovalores, Hopf teórico, barridos de a, b, tau e I |
| `src/visualization/plots.jl` | gráficos a partir de resultados ya calculados |
| `notebooks/01_individual.ipynb` | experimentación y visualización |
| `tests/runtests.jl` | validaciones analíticas y numéricas |
| `scripts/run_stage_A.jl` | generar figuras con un comando |
| `docs/ETAPA_A.md` | guía y resultados esperados |

La versión individual no implementa redes, controles ni ciclos límite por continuación. Eso corresponde a etapas posteriores.

### Nota sobre el caso b=0

La nulclina de la variable de recuperación deja de ser una función `w(v)` y pasa a ser la **recta vertical** `v=-a`. El código contempla expresamente este caso en el cálculo del equilibrio y los retratos de fase.
