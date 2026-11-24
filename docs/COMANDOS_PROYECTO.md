# Proyecto-COPTI — Guía de ejecución

**MA4703 · Control de sincronización en redes FitzHugh–Nagumo**

Esta guía reúne los comandos de las etapas A–G. Para la descripción matemática, los resultados, las figuras y sus limitaciones, consultar [RESULTADOS_GLOBALES.md](RESULTADOS_GLOBALES.md).

## 1. Preparar el entorno

Abrir la terminal de VS Code (PowerShell) en la raíz de `Proyecto-COPTI`, donde se encuentra `Project.toml`. Todos los comandos de esta guía se ejecutan desde esa carpeta.

```powershell
julia --version
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

Si se trabaja con una copia de GitHub y no hay cambios locales pendientes, actualizarla antes de ejecutar:

```powershell
git status
git pull --ff-only
```

## 2. Ejecutar las etapas

| Etapa | Contenido | Script | Pruebas |
|---|---|---|---|
| A | Dinámica individual, equilibrio y estabilidad | `run_stage_A.jl` | `runtests.jl` |
| B | Hopf, ciclos límite y shooting | `run_stage_B.jl` | `test_stage_B.jl` |
| C | Red, sincronización y estabilidad transversal | `run_stage_C.jl` | `test_stage_C.jl` |
| D | Controles suficientes y funciones de costo | `run_stage_D.jl` | `test_stage_D.jl` |
| E | Candidatos de control y diagnóstico de Pontryagin | `run_stage_E.jl` | `test_stage_E.jl` |
| F | Aproximaciones HJB y políticas de realimentación | `run_stage_F.jl` | `test_stage_F.jl` |
| G | Videos para la presentación | `run_stage_G_presentation_videos.jl` | — |

### A–F: simulaciones y experimentos

Ejecutar únicamente la línea de la etapa que se desea reproducir. Si se parte sin resultados previos, seguir el orden A–F.

```powershell
julia --project=. scripts/run_stage_A.jl
julia --project=. scripts/run_stage_B.jl
julia --project=. scripts/run_stage_C.jl
julia --project=. scripts/run_stage_D.jl
julia --project=. scripts/run_stage_E.jl
julia --project=. scripts/run_stage_F.jl
```

### Pruebas A–F

```powershell
julia --project=. tests/runtests.jl
julia --project=. tests/test_stage_B.jl
julia --project=. tests/test_stage_C.jl
julia --project=. tests/test_stage_D.jl
julia --project=. tests/test_stage_E.jl
julia --project=. tests/test_stage_F.jl
```

### G: videos de presentación

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl
```

Para regenerar los videos comparativos anteriores (basados en los candidatos de la etapa E), existe además:

```powershell
julia --project=. scripts/run_stage_F_video_only.jl
```

Si el script G indica que falta el paquete `Colors`, comprobar primero que se está usando el entorno del proyecto y ejecutar `Pkg.instantiate()`. Solo si sigue faltando como dependencia, añadirlo al entorno:

```powershell
julia --project=. -e 'using Pkg; Pkg.add("Colors")'
```

Esto puede modificar `Project.toml` y `Manifest.toml`.

## 3. Ejecutar los notebooks

Abrir la carpeta `notebooks/` en VS Code, seleccionar el kernel de Julia y ejecutar **Run All**. Si se usa Jupyter y falta el kernel de Julia, instalar IJulia en el entorno correspondiente (solo una vez):

```powershell
julia --project=. -e 'using Pkg; Pkg.add("IJulia")'
```

Para abrir Jupyter Notebook mediante IJulia (requiere Jupyter configurado):

```powershell
julia --project=. -e "using IJulia; notebook(dir=pwd())"
```

Si un notebook no encuentra el módulo o las dependencias, verificar que activa el proyecto raíz y reiniciar el kernel antes de volver a ejecutar las celdas.

Para listar los notebooks y scripts realmente presentes en la copia local:

```powershell
Get-ChildItem notebooks -Filter *.ipynb -Name
Get-ChildItem scripts -Filter *.jl -Name
```

## 4. Recopilar resultados

Con los experimentos A–F ya ejecutados, el recopilador consolida los resultados disponibles y ejecuta las pruebas, sin volver a realizar la optimización:

```powershell
julia --project=. scripts/recopilar_resultados.jl
```

## 5. Dónde encontrar las salidas

| Ruta | Contenido |
|---|---|
| `src/` | Módulos y funciones del modelo |
| `scripts/` | Scripts de ejecución y recopilación |
| `tests/` | Pruebas de las etapas A–F |
| `notebooks/` | Exploración y simulaciones interactivas |
| `data/processed/` | Resultados numéricos, principalmente JLD2 |
| `docs/figures/` | Figuras generadas |
| `docs/videos/` | Videos MP4 |
| `docs/RESULTADOS_GLOBALES.md` | Resultados, interpretación y conclusiones |

**Nota de mantenimiento:** los comandos A–F, pruebas, videos comparativos y recopilación se contrastaron con `RESULTADOS_GLOBALES.md`; el comando G corresponde al usado posteriormente para la presentación. No se verificó esta guía ejecutando los scripts ni contra una copia actual del repositorio. Si se renombró algún archivo, consultar los listados locales de `scripts/` y `notebooks/`.
