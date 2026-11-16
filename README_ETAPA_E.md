# MA4703 — Entrega Etapa E (se añade al repositorio A–D)

**No es otro repositorio.** Copia `src/`, `tests/`, `scripts/`, `notebooks/` y `docs/` sobre las carpetas existentes. Reemplaza solamente `src/FitzHughNagumoControl.jl`.

En Windows, desde la raíz del repositorio:

```powershell
julia --project=. tests/test_stage_E.jl
julia --project=. scripts/run_stage_E.jl
```

En `notebooks/05_pontryagin.ipynb` reinicia el kernel Julia y ejecuta las celdas en orden. No hay paquetes adicionales; se reutilizan las dependencias de A–D.

Esta entrega implementa **optimización directa aproximada y condiciones necesarias de Pontryagin** para ambos problemas. No se ha probado la optimalidad global y el control candidato puede no alcanzar la tolerancia terminal; se informa expresamente. Revisa [docs/ETAPA_E.md](docs/ETAPA_E.md) para fórmulas, algoritmos y limitaciones.
