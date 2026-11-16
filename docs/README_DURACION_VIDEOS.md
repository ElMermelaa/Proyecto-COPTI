# Etapa G — Corrección de duración

La Etapa E guarda 50 muestras del experimento de sincronización y aproximadamente 188 de desincronización. Un vídeo a 24 FPS con 50 cuadros dura solamente 2.08 segundos. La corrección no depende del número de muestras guardadas.

## Qué cambia

- `FPS_G = 12`
- `VIDEO_SECONDS_G = 20`
- `FINAL_HOLD_SECONDS_G = 4`
- Cada vídeo tendrá 240 fotogramas (aprox. 20 segundos).
- Los primeros 16 segundos recorren la trayectoria grabada mediante interpolación lineal entre estados consecutivos.
- Los últimos 4 segundos mantienen el último estado en pantalla para inspeccionar sus métricas.

**La pausa al final NO significa que se haya simulado el sistema más tiempo ni demuestra que ese estado se mantenga sin control.**

## Cómo instalar

Copia `scripts/run_stage_G_presentation_videos.jl` sobre el archivo homónimo en tu repositorio. No tienes que cambiar otros archivos ni instalar paquetes.

## Cómo ejecutar

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl --sync
julia --project=. scripts/run_stage_G_presentation_videos.jl --desync
```

O ambos:

```powershell
julia --project=. scripts/run_stage_G_presentation_videos.jl
```

Se regeneran los MP4 en `docs/videos/`.
