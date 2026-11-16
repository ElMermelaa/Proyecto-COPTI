# Etapa D · MA4703 · FitzHugh–Nagumo

Extensión incremental de las etapas A–C. Incluye controles suficientes por
acoplamiento y por corrientes individuales, certificados de Lyapunov,
corrientes diferenciadas para romper la igualdad, evaluación (no
optimización) de los costos J del cuaderno, pruebas, notebook, gráficos,
resultados en JLD2 y video opcional.

**IMPORTANTE:** copiar **el contenido** del ZIP en la raíz del repositorio.
Sólo se reemplaza `src/FitzHughNagumoControl.jl`, necesario para registrar
los módulos y exports de D. Los demás archivos son nuevos.

```powershell
julia --project=. tests/test_stage_D.jl
julia --project=. scripts/run_stage_D.jl
julia --project=. scripts/run_stage_D.jl --video
```

Consultar `docs/ETAPA_D.md` para las hipótesis matemáticas y los resultados.
