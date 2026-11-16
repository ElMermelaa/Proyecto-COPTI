# MA4703 — Etapa F: HJB y comparaciones

Esta entrega se **integra** en el repositorio de las etapas A–E; no crea uno nuevo.

1. Copia todo el **contenido** de este ZIP en `K:\Repositorios\Proyecto-COPTI`, combinando carpetas.
2. Reemplaza **solamente** `src/FitzHughNagumoControl.jl`. Los restantes archivos son nuevos.
3. Conserva `Project.toml`, `Manifest.toml`, tus `.jl` de A–E y resultados ya calculados.
4. Ejecuta `julia --project=. tests/test_stage_F.jl`.
5. Ejecuta `julia --project=. scripts/run_stage_F.jl`.
6. Para el video comparativo sin recalcular HJB: `julia --project=. scripts/run_stage_F_video_only.jl` (alternativamente, `run_stage_F.jl --video`).
7. Abre `notebooks/06_hjb.ipynb` y reinicia el kernel Julia.

Dependencias nuevas: **ninguna**. Reutiliza `CairoMakie`, `JLD2`, `OrdinaryDiffEq`, etc.

Hay dos aproximaciones HJB: sincronización exacta N=2 en subvariedad particular, y desincronización N=2 completa en grilla 4D, horizonte corto. **No deben interpretarse como soluciones óptimas globales del problema N=4 de la Etapa E.**

El archivo `docs/RESULTADOS_GLOBALES.md` contiene los resultados ya confirmados por tus ejecuciones hasta la Etapa E y los espacios por completar para F. **No borres imágenes previas**: se crean `F_*.png` y `F_*.mp4` nuevos.

No fue posible ejecutar Julia en el entorno donde se preparó este ZIP: antes de considerarlo validado necesitamos los resultados de `test_stage_F.jl` y del script en tu computador.
