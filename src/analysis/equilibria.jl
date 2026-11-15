"""    equilibrium(p; rtol=1e-13)

Devuelve `[v*, w*]`, el único equilibrio real del sistema.

- Si `b = 0`, `v* = -a` y `w* = f(-a)+I`.
- Si `b > 0`, se resuelve por bisección la ecuación estrictamente creciente

      (b/3)*v³ + (1-b)*v + a - b*I = 0.

Bisección aprovecha la unicidad demostrada en el cuaderno y evita elegir una
raíz espuria o un valor inicial delicado para Newton.
"""
function equilibrium(p::FHNParameters; rtol::Real=1e-13)
    (isfinite(rtol) && rtol > 0) ||
        throw(ArgumentError("rtol debe ser positivo y finito."))

    if iszero(p.b)
        vstar = -p.a
        return [vstar, fhn_nonlinearity(vstar) + p.I]
    end

    F(v) = (p.b / 3) * v^3 + (1 - p.b) * v + p.a - p.b * p.I

    left, right = -1.0, 1.0
    for _ in 1:1024
        F(left) <= 0 && break
        left *= 2
    end
    for _ in 1:1024
        F(right) >= 0 && break
        right *= 2
    end
    (F(left) <= 0 <= F(right)) ||
        error("No se pudo encerrar el equilibrio con aritmética Float64.")

    vstar = (left + right) / 2
    for _ in 1:250
        vstar = left + (right - left) / 2
        fmid = F(vstar)
        if fmid == 0 || (right - left) <= rtol * max(1.0, abs(vstar))
            break
        elseif fmid < 0
            left = vstar
        else
            right = vstar
        end
    end
    wstar = (vstar + p.a) / p.b
    return [vstar, wstar]
end

"""    equilibrium_residual(p, z=equilibrium(p))

Norma euclídea del campo vectorial en el candidato a equilibrio.
"""
equilibrium_residual(p::FHNParameters, z::AbstractVector=equilibrium(p)) =
    norm(fhn_rhs(z, p))
