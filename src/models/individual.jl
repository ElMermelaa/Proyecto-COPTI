"""    fhn_nonlinearity(v)

Función cúbica `f(v) = v - v³/3`.
"""
fhn_nonlinearity(v::Real) = v - v^3 / 3

"""    fhn_vector_field(v, w, p)

Campo vectorial `(dv/dt, dw/dt)` del modelo individual. No invierte tau.
"""
function fhn_vector_field(v::Real, w::Real, p::FHNParameters)
    dv = fhn_nonlinearity(v) - w + p.I
    dw = (v + p.a - p.b * w) / p.tau
    return (dv, dw)
end

"""    fhn_rhs(u, p)

Versión *out-of-place*, útil para diferencias automáticas con ForwardDiff.
"""
function fhn_rhs(u::AbstractVector, p::FHNParameters)
    length(u) == 2 || throw(ArgumentError("El estado debe tener dos componentes (v,w)."))
    dv, dw = fhn_vector_field(u[1], u[2], p)
    return [dv, dw]
end

"""    fhn!(du, u, p, t)

Versión *in-place* para `ODEProblem`, con `u = [v,w]`.
"""
function fhn!(du, u, p::FHNParameters, t)
    du[1], du[2] = fhn_vector_field(u[1], u[2], p)
    return nothing
end

"""    nullcline_v(v, p)

Nulclina de `dv/dt = 0`: `w = v-v³/3+I`.
"""
nullcline_v(v::Real, p::FHNParameters) = fhn_nonlinearity(v) + p.I

"""    nullcline_w(v, p)

Nulclina de `dw/dt = 0`, cuando `b > 0`: `w = (v+a)/b`.
Si `b = 0`, la nulclina es la recta VERTICAL `v = -a` y no
se puede escribir como función `w(v)`.
"""
function nullcline_w(v::Real, p::FHNParameters)
    iszero(p.b) && throw(DomainError(p.b,
        "Para b = 0, la nulclina de w es vertical: v = -a."))
    return (v + p.a) / p.b
end

"""    simulate_individual(p; u0=(-1.4,-0.2), tspan=(0.0,120.0), saveat=0.1,
                           solver=Tsit5(), reltol=1e-9, abstol=1e-9)

Integra la neurona individual y devuelve la `ODESolution` de SciML.
El estado se ordena `[v,w]`. Se usa un solucionador explícito por defecto;
si aparece rigidez, puede pasarse `solver=Rodas5P()`.
"""
function simulate_individual(p::FHNParameters;
        u0=(-1.4, -0.2),
        tspan=(0.0, 120.0),
        saveat=0.1,
        solver=Tsit5(),
        reltol::Real=1e-9,
        abstol::Real=1e-9)

    length(u0) == 2 || throw(ArgumentError("u0 debe tener dos componentes (v,w)."))
    all(x -> x isa Real && isfinite(x), u0) ||
        throw(ArgumentError("u0 debe contener dos números reales finitos."))
    length(tspan) == 2 || throw(ArgumentError("tspan debe ser (t0, tf)."))
    t0, tf = Float64.(tspan)
    (isfinite(t0) && isfinite(tf) && tf > t0) ||
        throw(ArgumentError("Se exige t0 < tf y tiempos finitos."))
    if saveat isa Real
        (isfinite(saveat) && saveat > 0) ||
            throw(ArgumentError("saveat debe ser positivo y finito."))
    end
    (isfinite(reltol) && reltol > 0 && isfinite(abstol) && abstol > 0) ||
        throw(ArgumentError("Las tolerancias deben ser positivas y finitas."))

    prob = ODEProblem(fhn!, Float64.(collect(u0)), (t0, tf), p)
    return solve(prob, solver; saveat=saveat, reltol=reltol, abstol=abstol)
end

"""    simulate_ensemble(p, initial_conditions; kwargs...)

Integra una colección de condiciones iniciales sin mezclarlas en un único modelo.
Devuelve un vector de soluciones, útil para retratos de fase.
"""
function simulate_ensemble(p::FHNParameters, initial_conditions; kwargs...)
    return [simulate_individual(p; u0=u0, kwargs...) for u0 in initial_conditions]
end
