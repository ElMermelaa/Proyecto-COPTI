# MA4703 — Resumen automático de resultados A–F

Generado: 2026-11-15 20:32:37

El reporte recoge archivos existentes. No demuestra optimalidad de PMP ni convergencia de HJB.

## Validación

| Etapa | Resultado | Tiempo (s) | Log |
|---|---|---:|---|
| A | Aprobada | 13.8 | tests/test_A.txt |
| B | Aprobada | 14.9 | tests/test_B.txt |
| C | Aprobada | 27.8 | tests/test_C.txt |
| D | Aprobada | 28.2 | tests/test_D.txt |
| E | Aprobada | 39.3 | tests/test_E.txt |
| F | Aprobada | 12.0 | tests/test_F.txt |

## Archivos JLD2

| Archivo | Tamaño | Campos | Errores de lectura |
|---|---:|---:|---:|
| data/processed/C_notebook_lambda_m1p5.jld2 | 236.2 KB | 21 | 0 |
| data/processed/etapa_B_resultados.jld2 | 23.3 KB | 17 | 0 |
| data/processed/etapa_C_lambda_0p45.jld2 | 236.3 KB | 21 | 0 |
| data/processed/etapa_C_lambda_1p5.jld2 | 236.3 KB | 21 | 0 |
| data/processed/etapa_C_lambda_m0p3.jld2 | 236.3 KB | 21 | 0 |
| data/processed/etapa_C_lambda_m1p5.jld2 | 236.3 KB | 21 | 0 |
| data/processed/etapa_D_desync_0p0.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_desync_0p25.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_desync_0p75.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_feedback.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_sync_0p45.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_sync_m0p15.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_sync_m0p3.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_sync_m1p5.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_D_sync_modulada.jld2 | 334.4 KB | 27 | 0 |
| data/processed/etapa_E_desync.jld2 | 32.9 KB | 23 | 0 |
| data/processed/etapa_E_sync.jld2 | 15.4 KB | 23 | 0 |
| data/processed/etapa_F_hjb.jld2 | 75.0 KB | 28 | 0 |

Todos los campos y rangos de datos están en `datos/`. Las series temporales reducidas están en `curvas/`.

## Gráficos y videos

- Figuras PNG encontradas: **29**; incluidas en el paquete: **29**.

- Videos/GIF encontrados: **5**. No se adjuntan automáticamente para evitar un ZIP gigantesco.

- Archivos de inventario: `inventario_figuras.csv` e `inventario_videos.csv`.

## Precauciones de interpretación

- La condición `E(T) ≤ δ` debe verificarse independientemente en cada candidato de sincronización.

- Un residuo PMP grande impide presentar un candidato como extremal verificado.

- Una HJB fuera de cobertura no suministra una política numérica validada.

- Solo comparar costos cuando `N`, parámetros, condiciones iniciales, restricciones y pesos son iguales.
