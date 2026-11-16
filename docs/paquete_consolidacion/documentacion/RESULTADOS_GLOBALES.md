# MA4703 — Registro global de resultados

> Borrador consolidado. Completar con los valores de la Etapa F una vez que el script se haya ejecutado en Julia. Los costos de experimentos con distinto N, estados iniciales, horizontes o pesos **no son directamente comparables**.

## Estado de implementación

| Etapa | Qué se implementó | Pruebas reportadas | Evidencia e interpretación |
|---|---|---:|---|
| A | Modelo individual, equilibrio y estabilidad, nulclinas y gráficos | 61/61 | Para `a=.7,b=.5,τ=12.5,I=1.4`, equilibrio `(0,1.4)` inestable |
| B | Hopf, ciclos límite, shooting y Floquet | 38/38 | `I₋≈0.10667`, `I₊≈2.69333`; ciclo `T≈34.5696` en `I=1.4` |
| C | Red N=12, E(t), r(t), gráficos y video neuronal | 89/89 | `λ<0` atractivo; `λ>0` puede separar estados |
| D | Leyes suficientes, Lyapunov, evaluaciones de J y corrientes diferenciadas | 161/161 | Se observó sincronía para λ negativo incluso fuera del criterio suficiente |
| E | Controles candidatos, condiciones de PMP y diagnósticos | 54/54 | Sincronización factible, PMP todavía no verificado con residuo pequeño |
| F | HJB estacionaria reducida N=2, temporal N=2 y comparaciones | **Pendiente de ejecución** | Grillas y controles exploratorios; no extrapolar a N=4 |

## Comparación E: sincronización (N=4, mismo estado inicial)

| Estrategia | Primera llegada a E≤0.03 | J de llegada | Restricción |
|---|---:|---:|---|
| λ=-1.5 constante | 10.29036 | 21.66269 | Cumple |
| λ(t) candidata Etapa E | 7.25873 | 14.61266 | Cumple |

**Advertencias:** El control variable es factible, pero los residuos informados de PMP fueron aprox. `0.9984` (ley proyectada) y `1.1922` (Hamiltoniano terminal), por lo que **no** se certifica como extremal de Pontryagin ni como óptimo global. La desigualdad terminal real verificada con Tsit5 fue `E(T)≈0.02989966`, con margen numérico respecto de δ.

## Comparación E: desincronización (N=4, horizonte fijo)

| Estrategia | J | r²(T) | Observación |
|---|---:|---:|---|
| u≡0 | 92.0000 | 1 si permanece sincronizada | Línea base |
| Candidato de corrientes | 30.2381 (Tsit5) | 0.00063582 | Suma de corrientes ≈0; amplitudes ≤0.8 |

**Advertencias:** `r²≈0` no prueba ausencia de dos clústeres opuestos ni desincronización sostenida. Residuo PMP de corrientes reportado ≈`1.0516`; no se certifica optimalidad.

## Resultados Etapa F (completar tras ejecutar)

| Experimento | J HJB | J referencia | Cobertura de malla | Nota |
|---|---:|---:|---|---|
| F1 llegada sincronización exacta reducida N=2 | Pendiente | Pendiente (λ=-1.5) | Pendiente | `r²=0` fuera del objetivo por geometría antipodal |
| F2 desincronización en 4D N=2 | Pendiente | Pendiente (`u=0`) | Pendiente | Horizonte corto, acciones discretizadas |

## Recursos y archivos

- `docs/figures/A_*.png` a `F_*.png`: figuras por etapa (si existen).
- `docs/videos/C_*.mp4`, `D_*.mp4`: experimentos visuales previos.
- `docs/videos/F_01_comparacion_control_variable_vs_constante.mp4`: video opcional desde sincronización E.
- `docs/videos/F_02_comparacion_corrientes_vs_sin_control.mp4`: video opcional desde desincronización E.
- `data/processed/etapa_B_resultados.jld2`, `etapa_E_sync.jld2`, `etapa_E_desync.jld2`, `etapa_F_hjb.jld2`.

### Pendientes para el informe final

1. Validar y anotar resultados reales del script F.
2. Estudiar convergencia de grilla para HJB.
3. Refinar extremales de Pontryagin si hay tiempo.
4. Comparar controles **solamente** sobre el mismo modelo, pesos, N, horizonte y estado inicial.
5. Interpretar por separado sincronización de estados, alineamiento de fase y persistencia de oscilaciones.
