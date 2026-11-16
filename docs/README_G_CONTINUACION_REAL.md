# Corrección definitiva: continuar la dinámica después del control óptimo candidato

Esta versión **sí prolonga la simulación física**, en vez de ralentizar o congelar la reproducción.

## Instalación

Reemplaza **solo** `scripts/run_stage_G_presentation_videos.jl` en el repositorio. No requiere paquetes nuevos ni repetir optimizaciones.

## Ejecutar

Desde la raíz del repositorio:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl
```

O un solo video:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl --sync
julia --project=. scripts/run_stage_G_presentation_videos.jl --desync
```

Los videos nuevos se guardan como:

- `docs/videos/G_01_presentacion_sincronizacion_extendido.mp4`
- `docs/videos/G_02_presentacion_desincronizacion_extendido.mp4`

**No sobrescribe los videos anteriores**.

## Qué pasa después de T

Las políticas de Etapa E solo fueron calculadas hasta el tiempo final T. Por ello se debe especificar una ley posterior. Por defecto, el programa **mantiene constante el último valor del control candidato**, tanto para lambda como para las corrientes individuales. Esto es una continuación controlada, **no una afirmación de que ese control sea óptimo después de T**.

Si prefieres examinar qué ocurre al **retirar el control** después de T:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl --release
```

Esto pone lambda=0 en sincronización y u_i=0 en desincronización; conserva lambda0=-0.15 en este último caso. Se generan videos distintos con sufijo `_libre.mp4`.

## Cómo aumentar el tiempo físico

En las primeras líneas del script:

```julia
const POST_SYNC_G = 80.0      # unidades temporales físicas después de T
const POST_DESYNC_G = 80.0    # unidades temporales físicas después de T
const FPS_G = 24             # FPS de reproducción
const VIDEO_SECONDS_G = 24   # duración de MP4
```

Por ejemplo, `POST_SYNC_G=120.0` simulará 120 unidades temporales **adicionales** desde el estado final de la Etapa E. El video recorrerá esa trayectoria física real, no la reproducirá más lentamente.

Se marca T con una línea roja discontinua en los gráficos, y el pie distingue cuándo se pasa a continuación con control constante o retirado.

## Advertencia científica

Que las neuronas se mantengan sincronizadas o desincronizadas después de T **depende de la ley posterior** y no se garantiza por la optimización anterior. La persistencia observada es un resultado numérico exploratorio.
