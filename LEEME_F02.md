# Correccion diagnostica F_02 (Etapa F)

## Que problema corrige

La imagen `docs/figures/F_02_hjb_desync_corte_4D.png` podia salir vacia porque las celdas del corte de la funcion de valor en tiempo 0 eran `Inf` (falta de cobertura de la malla 4D). Ademas, el codigo anterior seleccionaba por accidente la primera accion de la lista (normalmente `-Umax`) cuando TODAS las acciones eran invalidas. Esa accion NO se deducia de Bellman.

Esta correccion evita ese control arbitrario y sustituye, **solo para diagnostico**, por la corriente neutra `u=(0,0)`, siempre con `covered=false`. La nueva figura presenta el mapa de cobertura o el valor disponible, el costo terminal, `r²(t)` y la corriente efectivamente aplicada. Cuando HJB carece de cobertura, lo indica explicitamente. **No resuelve ni valida la HJB de cuatro variables.**

## Instalacion

1. Extraer ZIP. Copiar *el contenido* de las carpetas `src/` y `scripts/` en la raiz de `K:\Repositorios\Proyecto-COPTI`.
2. Aceptar fusion de carpetas; reemplazar **solo**:
   - `src/control/hjb_desync_small.jl`
   - `src/visualization/hjb_plots.jl`
   - `scripts/run_stage_F.jl`
3. **No borrar** imágenes ni datos A-E; la imagen `F_02_hjb_desync_corte_4D.png` se regenerara automaticamente.
4. En terminal del repositorio:
   ```powershell
   julia --project=. tests/test_stage_F.jl
   julia --project=. scripts/run_stage_F.jl
   ```

No hay nuevas dependencias. Esta version NO fue ejecutada con Julia en el entorno de preparacion: las pruebas deben correrse localmente.

## Interpretacion

Si el nuevo grafico informa cobertura insuficiente, eso sigue siendo una limitacion real del experimento: se necesita mejorar malla, horizonte, condiciones de borde o el algoritmo antes de concluir sobre realimentacion HJB. El caso N=2,T=1.2 tampoco muestra una disminucion apreciable de `r²` con los parametros actuales. No presentar como control HJB optimo.

Los videos `F_01_*.mp4` y `F_02_*.mp4` comparan controles de Etapa E, no HJB contra Pontryagin; esta correccion no los modifica.
