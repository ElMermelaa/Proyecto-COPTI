# MA4703 — Redes FitzHugh–Nagumo y control óptimo
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