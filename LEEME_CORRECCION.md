# Corrección Etapa E — restricción terminal E(T) ≤ δ

Se verificó en el cuaderno teórico que la restricción **siempre estuvo en la formulación** de sincronización. El primer optimizador la imponía mediante un Lagrangiano aumentado, pero aceptaba candidatos infactibles cuando se agotaban las iteraciones. Esta corrección **no cambia el problema matemático**.

## Archivos que hay que reemplazar

1. `src/control/optimal_direct.jl` — mantiene el Lagrangiano aumentado, pero añade recuperación explícita de factibilidad: controles constantes testigo, búsqueda de llegada, interpolación y refinamiento sólo con candidatos factibles. Comprueba RK4 y Tsit5. Si no encuentra solución factible, lo reporta claramente.
2. `scripts/run_stage_E.jl` — detiene la creación de resultados de sincronización si no se cumple la restricción.
3. `tests/test_stage_E.jl` — comprueba la coherencia del indicador `feasible` y añade una prueba de factibilidad.
4. `notebooks/05_pontryagin.ipynb` — mismas comprobaciones para Jupyter.

**No reemplazar** `src/FitzHughNagumoControl.jl`, archivos A–D ni `Project.toml`/`Manifest.toml`.

## Ejecutar

```powershell
julia --project=. tests/test_stage_E.jl
julia --project=. scripts/run_stage_E.jl
```

Reinicia el kernel antes de volver a abrir el notebook actualizado.

## Qué borrar

**No es necesario borrar** figuras ni datos: `run_stage_E.jl` sobrescribe `docs/figures/E_01_sincronizacion_candidata.png`, `E_02_residuos_PMP_sync.png` y `data/processed/etapa_E_sync.jld2` cuando la verificación sale bien. Los archivos antiguos de A–D se conservan; tampoco es necesario eliminar los de desincronización. Hasta que la nueva ejecución sea satisfactoria, **no uses** las figuras antiguas de sincronización E para justificar conclusiones.

## Advertencias

- El método sigue siendo una aproximación de control constante a trozos. Recuperar factibilidad no prueba que la solución sea un mínimo local, ni que satisfaga con precisión Pontryagin.
- El multiplicador `mu` es aproximado, y los residuos PMP deben revisarse de nuevo tras la corrección.
- Si con las cotas y el horizonte elegido no existe un control factible, `feasible=false`: el código **no inventa** una solución.
- La condición terminal de un óptimo libre suele ser activa, `E(T)=δ`. La búsqueda numérica usa un pequeño margen interior para evitar que el error de integración viole la desigualdad.
