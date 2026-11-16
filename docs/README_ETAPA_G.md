# Etapa G — Videos de presentación (versión corregida)

La versión anterior dependía de `Colors` y `ColorSchemes` sin registrarlos en el entorno y no reconocía los nombres de campos de la Etapa E. Esta versión utiliza solamente **CairoMakie y JLD2**, ya disponibles en el repositorio. 

No modifica resultados anteriores ni vuelve a optimizar.

## Instalación

Reemplazar **únicamente** `scripts/run_stage_G_presentation_videos.jl`.

## Ejecución

En la raíz del repositorio:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl
```

Puedes generar un solo video:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl --sync
julia --project=. scripts/run_stage_G_presentation_videos.jl --desync
```

## Entrada

- `data/processed/etapa_E_sync.jld2`
- `data/processed/etapa_E_desync.jld2`

## Salida

- `docs/videos/G_01_presentacion_sincronizacion.mp4`
- `docs/videos/G_02_presentacion_desincronizacion.mp4`

En cada video se presentan: red neuronal pulsante, plano de fase `(v,w)`, circunferencia de fases geométricas, `E(t)`, `r(t)` y el control aplicado. Todos se basan en **los mismos datos guardados de la Etapa E**.

**Precaución:** los controles de la Etapa E son candidatos numéricos factibles, no óptimos globales certificados. Si la distancia de una neurona al centro geométrico se hace cero, su fase no está definida.

No se ha ejecutado Julia en el entorno de preparación de este paquete: la comprobación funcional final debe realizarse en tu PC.
