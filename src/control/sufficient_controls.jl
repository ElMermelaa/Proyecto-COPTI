# ETAPA D — Controles suficientes y perturbaciones diferenciadas.
# Depende de FHNNetwork, sin alterar ni reimplementar el campo de la Etapa C.

"""    AttractiveModulatedCoupling(ceiling, depth, omega)

Familia de controles temporales λ(t)=ceiling-depth*sin(omega*t)^2.
Como depth>=0, el SUPREMO GLOBAL es <= ceiling (conocido analíticamente).
Esto permite un certificado para λ(t) no constante sin inferir cotas desde
muestras discretas. Si ceiling<-1 y b>0 hay sincronía garantizada.
"""
struct AttractiveModulatedCoupling
    ceiling::Float64
    depth::Float64
    omega::Float64
    function AttractiveModulatedCoupling(ceiling::Real,depth::Real,omega::Real)
        all(isfinite,(ceiling,depth,omega)) && depth>=0 && omega>=0 ||
            throw(ArgumentError("Parámetros no finitos o depth/omega negativos."))
        new(Float64(ceiling),Float64(depth),Float64(omega))
    end
end

(c::AttractiveModulatedCoupling)(t::Real) =
    c.ceiling - c.depth*sin(c.omega*t)^2

"""    feedback_currents(x, k)

Realimentación del cuaderno: `u_i = -k (v_i - v̄)`, aplicada a la
ecuación del potencial. Siempre tiene suma cero. No aplica saturaciones.
"""
function feedback_currents(x::AbstractVector, k::Real)
    N = _validate_state_size(x)
    isfinite(k) && k >= 0 ||
        throw(ArgumentError("La ganancia k debe ser finita y no negativa."))
    vm = sum(x[1:2:end]) / N
    return [-k*(x[2i-1]-vm) for i in 1:N]
end

"""    FHNCurrentFeedback(base::FHNNetwork; k=1.5)

Sistema de la Etapa C MÁS `u_i(x) = -k(v_i-v̄)`.
`base.currents` representa las corrientes adicionales preexistentes; el
feedback se suma a ellas. Para las garantías de Lyapunov, deben ser comunes.
"""
struct FHNCurrentFeedback{NET<:FHNNetwork}
    base::NET
    k::Float64
end

function FHNCurrentFeedback(base::FHNNetwork; k::Real=1.5)
    isfinite(k) && k >= 0 ||
        throw(ArgumentError("k debe ser finita y >= 0."))
    return FHNCurrentFeedback{typeof(base)}(base, Float64(k))
end

"""    feedback_rhs!(du, x, system, t)

Calcula la dinámica EXACTA bajo la realimentación sin modificar `network_rhs!`.
"""
function feedback_rhs!(du, x, system::FHNCurrentFeedback, t)
    network_rhs!(du, x, system.base, t)
    for (i, ui) in enumerate(feedback_currents(x, system.k))
        du[2i-1] += ui
    end
    return nothing
end

function feedback_rhs(x::AbstractVector, system::FHNCurrentFeedback, t::Real=0.0)
    out = similar(x)
    feedback_rhs!(out, x, system, t)
    return out
end

"""    simulate_feedback(system; u0, tspan, saveat, solver, reltol, abstol)

Integra la red bajo corrientes dependientes del estado. Devuelve ODESolution.
"""
function simulate_feedback(system::FHNCurrentFeedback;
        u0=distributed_initial_state(system.base.N, system.base.neuron),
        tspan=(0.0, 200.0), saveat=0.2, solver=Tsit5(),
        reltol::Real=1e-9, abstol::Real=1e-9)
    length(u0) == 2 * system.base.N ||
        throw(ArgumentError("El estado inicial debe tener 2N componentes."))
    all(q -> q isa Real && isfinite(q), u0) ||
        throw(ArgumentError("Estado inicial inválido."))
    length(tspan) == 2 || throw(ArgumentError("tspan debe tener 2 entradas."))
    t0, tf = Float64.(tspan)
    all(isfinite, (t0,tf)) && tf > t0 ||
        throw(ArgumentError("El intervalo de tiempo debe ser creciente."))
    saveat isa Real && isfinite(saveat) && saveat > 0 ||
        throw(ArgumentError("saveat debe ser positivo."))
    isfinite(reltol) && reltol > 0 && isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("Tolerancias no válidas."))
    prob = ODEProblem(feedback_rhs!, Float64.(collect(u0)), (t0,tf), system)
    return solve(prob, solver; saveat=saveat, reltol=reltol, abstol=abstol)
end

"""    balanced_current_profile(N; pattern=:split)

Perfil determinista con suma cero y norma infinita 1. `:split` divide la red
en dos grupos con estímulos de signo contrario; `:alternating` alterna nodos;
`:sinusoidal` crea un patrón suave. No garantiza desincronización sostenida.
"""
function balanced_current_profile(N::Integer; pattern::Symbol=:split)
    N >= 2 || throw(ArgumentError("Se requieren al menos dos neuronas."))
    raw = if pattern == :split
        [i <= fld(N,2) ? 1.0 : -1.0 for i in 1:N]
    elseif pattern == :alternating
        [isodd(i) ? 1.0 : -1.0 for i in 1:N]
    elseif pattern == :sinusoidal
        [sin(2pi*(i-1)/N) for i in 1:N]
    else
        throw(ArgumentError("Patrón desconocido: $pattern"))
    end
    raw .-= sum(raw)/N
    m = maximum(abs, raw)
    m > 1e-12 || throw(ArgumentError("Patrón degenerado."))
    raw ./= m
    # Compensación del redondeo: suma prácticamente cero incluso con N impar.
    raw .-= sum(raw)/N
    return raw
end

"""    balanced_current_law(profile; amplitude=0.5, ramp_time=0.0)

Devuelve `(t -> u(t))`, con suma_i u_i(t) ≈ 0 y |u_i(t)| <= amplitude.
Con ramp_time=0 el estímulo es constante desde t=0; si es positivo,
se enciende suavemente como 1-exp(-max(t,0)/ramp_time).
"""
function balanced_current_law(profile::AbstractVector;
        amplitude::Real=0.5, ramp_time::Real=0.0)
    N = length(profile)
    N >= 2 && all(q -> q isa Real && isfinite(q), profile) ||
        throw(ArgumentError("Perfil real, finito y de longitud >= 2 requerido."))
    abs(sum(profile)) <= 1e-10 ||
        throw(ArgumentError("El perfil debe tener suma cero."))
    amplitude >= 0 && isfinite(amplitude) ||
        throw(ArgumentError("amplitude debe ser finita y no negativa."))
    ramp_time >= 0 && isfinite(ramp_time) ||
        throw(ArgumentError("ramp_time debe ser finito y no negativo."))
    m = maximum(abs, profile)
    shape = m == 0 ? zeros(Float64, N) : Float64.(profile)./m
    return t -> begin
        g = ramp_time == 0 ? 1.0 : -expm1(-max(Float64(t),0.0)/ramp_time)
        amplitude * g .* shape
    end
end

"""    control_current_matrix(sol, sys)

Matriz N×nt de corrientes INDIVIDUALES adicionales a I0 (incluye feedback).
No incluye la corriente basal `p.I`.
"""
function control_current_matrix(sol, net::FHNNetwork)
    N = net.N
    out = zeros(Float64, N, length(sol.t))
    for (j,t) in enumerate(sol.t)
        u = _current_at(net.currents, t)
        length(u) == N || throw(DimensionMismatch("Corrientes de dimensión incorrecta."))
        out[:,j] .= u
    end
    return out
end

function control_current_matrix(sol, system::FHNCurrentFeedback)
    U = control_current_matrix(sol, system.base)
    for (j,x) in enumerate(sol.u)
        U[:,j] .+= feedback_currents(x, system.k)
    end
    return U
end

"""    total_current_effort(sol, sys)

Integral numérica ∫sum_i u_i(t)^2 dt de las corrientes adicionales.
Esta métrica es descriptiva; NO reemplaza el costo de sincronización en λ.
"""
function total_current_effort(sol, sys)
    U = control_current_matrix(sol, sys)
    integrand = vec(sum(abs2, U; dims=1))
    return _trapezoidal_integral(sol.t, integrand)
end
