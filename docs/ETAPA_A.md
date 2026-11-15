# Etapa A — Dinámica individual: fundamentos, validación y experimentos

## 1. Ecuaciones oficiales

\[
\dot v = f(v)-w+I,\qquad \dot w=(v+a-bw)/\tau,
\quad f(v)=v-v^3/3.
\]

Hipótesis del enunciado: `a ≥ 0`, `0 ≤ b < 1`, `tau > 0`, corriente externa real `I`.

**Esta etapa no incluye acoplamiento**; el signo de `+lambda (vi-v̄)` será importante al implementar las redes, pero aquí estudiamos una neurona aislada.

## 2. Equilibrio único (cuaderno, capítulo 1)

En un equilibrio tenemos

\[
0=f(v^*)-w^*+I,\qquad 0=v^*+a-bw^*.
\]

- Para `b > 0`: `w* = (v*+a)/b` y `v*` es la única raíz real de `(b/3)v³ + (1-b)v + a-bI = 0`. Su derivada es `b v²+(1-b)>0`, así que usamos bisección con encierro progresivo.
- Para `b = 0`: `v*=-a`, `w*=f(-a)+I`. No se divide por `b`.

La función `equilibrium_residual` evalúa `norm(F(v*,w*))` para comprobar la solución.

## 3. Nulclinas y jacobiano

- `dv/dt=0` equivale a `w=f(v)+I` (cúbica).
- `dw/dt=0` equivale a `w=(v+a)/b` si `b>0`, y a `v=-a` si `b=0`.

\[
J(v,w)=\begin{pmatrix} 1-v^2 & -1\\ 1/\tau & -b/\tau\end{pmatrix},
\quad \operatorname{tr}J^*=1-(v^*)^2-\frac b\tau,
\quad \det J^*=\frac{1-b+b(v^*)^2}{\tau}>0.
\]

Por el criterio traza-determinante: traza negativa = equilibrio localmente asintóticamente estable; traza positiva = inestable; traza nula = la linealización no decide estabilidad no lineal. `stability_info` devuelve eigenvalores y clasificación.

## 4. Currientes candidatas de Hopf (cuaderno, capítulo 2)

Solo si `0 < b < 1` y `tau > b`:

\[
s=\sqrt{1-b/\tau},\qquad
v_H^\pm=\pm s,\qquad
I_H^\pm=\frac{(v_H^\pm)^3}{3}-v_H^\pm+\frac{v_H^\pm+a}{b}.
\]

`hopf_currents` calcula estos valores y `omega = sqrt(tau-b²)/tau`. Esta etapa marca el cambio de signo de la traza; **no** sigue ramas de ciclos ni prueba numéricamente su estabilidad. El estudio de bifurcaciones y de órbitas periódicas es la Etapa B.

## 5. Ejemplo numérico sugerido

Valores elegidos para demostrar la API, no para alterar el modelo teórico:

```julia
p = FHNParameters(a=0.7, b=0.5, tau=12.5, I=1.4)
```

Se tiene exactamente `v*=0`, `w*=1.4`, `tr(J*)=0.96>0`, por lo que el equilibrio es inestable. El modelo tiene un régimen donde se observan oscilaciones para muchas condiciones iniciales, pero una gráfica aislada **no certifica** un ciclo límite. El cálculo de su período y atractividad es trabajo para la Etapa B.

Para explorar equilibrios estables, variar `I` fuera del intervalo crítico reportado por `hopf_currents(p)`.

## 6. Qué comprobar en las figuras

1. **Series temporales:** `v(t)` y `w(t)` desde una condición inicial, mostrando oscilaciones, transitorios o aproximación a un equilibrio según el caso.
2. **Retrato de fase:** cúbica roja (`dv/dt=0`), recta verde (`dw/dt=0`, vertical si `b=0`), punto fijo negro, campo de direcciones y varias trayectorias.
3. **Barrido de corriente:** curva `v*(I)`, traza de `J*` y marcas de corrientes candidatas de Hopf.
4. **Validación jacobiano:** el jacobiano analítico debe coincidir con `ForwardDiff.jacobian` dentro de tolerancia numérica.
5. **Variación paramétrica:** observar cómo cambian el equilibrio y la traza al variar individualmente `a`, `b`, `tau` e `I`, sin suponer que la estabilidad local describe toda la dinámica.
6. **Validación integración:** una trayectoria que parte exactamente del equilibrio permanece allí y dos integraciones con distintas tolerancias deben coincidir aproximadamente.

## 7. Qué NO se afirma aún

- No se deduce globalidad de un ciclo estable de ver una trayectoria cerrada.
- No se demuestra existencia o tipo de bifurcación solo por `tr(J)=0`.
- No se usan automáticamente resultados de Hopf cuando `b=0` o `tau≤b`.
- No se confunde estabilidad de un equilibrio individual con estabilidad de la sincronización de la red.

## 8. Comandos desde la raíz del repositorio

```powershell
julia --project=. tests/runtests.jl
julia --project=. scripts/run_stage_A.jl
```

Si al cargar un notebook ya ejecutado editas `src/`, reinicia el kernel de Julia y vuelve a ejecutar las celdas para utilizar el código actualizado.
