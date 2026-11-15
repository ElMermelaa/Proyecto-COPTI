# Métricas de dispersión y fase del cuaderno (caps. 3–4).

"""    synchronization_error(x)

E(x) = (1/N) * Σᵢ ((vᵢ-v̄)² + (wᵢ-w̄)²).
Se anula si y solo si todos los estados son idénticos.
"""
function synchronization_error(x::AbstractVector)
    N = _validate_state_size(x)
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    return sum((x[2i-1]-vm)^2 + (x[2i]-wm)^2 for i in 1:N)/N
end

"""    lyapunov_dispersion(x, tau)

V(x) = (1/2) Σᵢ (vᵢ-v̄)² + (τ/2) Σᵢ (wᵢ-w̄)².
"""
function lyapunov_dispersion(x::AbstractVector, tau::Real)
    N = _validate_state_size(x)
    isfinite(tau) && tau > 0 || throw(ArgumentError("tau debe ser positivo."))
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    return 0.5*sum((x[2i-1]-vm)^2 + tau*(x[2i]-wm)^2 for i in 1:N)
end

function _validate_state_size(x)
    iseven(length(x)) && length(x) >= 4 ||
        throw(DimensionMismatch("Se requieren 2N componentes, N>=2."))
    all(z -> z isa Real && isfinite(z), x) ||
        throw(ArgumentError("El estado debe ser real y finito."))
    return div(length(x), 2)
end

"""    phase_order_squared(x, center; center_tol=1e-10, strict=false)

qᵢ = (zᵢ-c)/||zᵢ-c||; r²=||(1/N)Σᵢqᵢ||².
Si alguna neurona se acerca a `center` a distancia <= `center_tol`, la fase
NO está definida: se devuelve NaN, o se lanza DomainError si strict=true.
No se fuerza r=0 o r=1 arbitrariamente.
"""
function phase_order_squared(x::AbstractVector, center;
        center_tol::Real=1e-10, strict::Bool=false)
    N = _validate_state_size(x)
    length(center) == 2 && all(isfinite, center) ||
        throw(ArgumentError("center debe contener dos números finitos."))
    isfinite(center_tol) && center_tol >= 0 ||
        throw(ArgumentError("center_tol debe ser finito y no negativo."))
    qv = 0.0
    qw = 0.0
    for i in 1:N
        dv = x[2i-1] - center[1]
        dw = x[2i] - center[2]
        dist = hypot(dv, dw)
        if dist <= center_tol
            strict && throw(DomainError(dist, "Fase indefinida: neurona cerca del centro."))
            return NaN
        end
        qv += dv / dist
        qw += dw / dist
    end
    return clamp((qv*qv + qw*qw)/N^2, 0.0, 1.0)
end

"""    phase_order(x, center; kwargs...)

Devuelve r = sqrt(r²), o NaN cuando la fase no está definida.
"""
phase_order(x::AbstractVector, center; kwargs...) =
    sqrt(phase_order_squared(x, center; kwargs...))

"""    network_metrics(sol, net; center=equilibrium(net.neuron), center_tol=1e-10)

Trazas E(t), V(t), r²(t), r(t), distancia mínima al centro, y máscara
`phase_defined`. No normaliza trayectorias ni modifica la solución.
"""
function network_metrics(sol, net::FHNNetwork;
        center=equilibrium(net.neuron), center_tol::Real=1e-10)
    N = net.N
    t = Float64.(sol.t)
    E = Float64[]
    V = Float64[]
    r2 = Float64[]
    mind = Float64[]
    for state in sol.u
        length(state) == 2N || throw(DimensionMismatch("Solución de dimensión incorrecta."))
        push!(E, synchronization_error(state))
        push!(V, lyapunov_dispersion(state, net.neuron.tau))
        push!(r2, phase_order_squared(state, center; center_tol=center_tol))
        push!(mind, minimum(hypot(state[2i-1]-center[1],
                                 state[2i]-center[2]) for i in 1:N))
    end
    return (t=t, E=E, V=V, r2=r2, r=sqrt.(r2),
            min_center_distance=mind, phase_defined=isfinite.(r2),
            center=Float64.(collect(center)))
end

"""    is_mean_field(net; atol=1e-12)

Verifica si Aᵢⱼ=1 fuera de la diagonal (la diagonal no afecta la dinámica).
"""
function is_mean_field(net::FHNNetwork; atol::Real=1e-12)
    isfinite(atol) && atol >= 0 || throw(ArgumentError("atol debe ser no negativo."))
    for i in 1:net.N, j in 1:net.N
        if i != j && abs(net.A[i,j] - 1.0) > atol
            return false
        end
    end
    return true
end

"""    lyapunov_derivative(x, net, t=0.0)

Derivada exacta de V calculada con el campo de la EDO, no por diferencia finita.
"""
function lyapunov_derivative(x::AbstractVector, net::FHNNetwork, t::Real=0.0)
    length(x) == 2 * net.N || throw(DimensionMismatch("Estado incorrecto."))
    dx = network_rhs(x, net, t)
    N = net.N
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    return sum((x[2i-1]-vm)*dx[2i-1] +
               net.neuron.tau*(x[2i]-wm)*dx[2i] for i in 1:N)
end

"""    lyapunov_upper_bound(x, net, t=0.0)

Lado derecho de V̇ <= (1+λ)Σeᵢ² - bΣhᵢ², VÁLIDO solamente para
A todos-con-todos y corriente individual común. En otro caso, error.
"""
function lyapunov_upper_bound(x::AbstractVector, net::FHNNetwork, t::Real=0.0)
    is_mean_field(net) ||
        throw(ArgumentError("La cota del cuaderno requiere campo medio completo."))
    u = _current_at(net.currents, t)
    u isa AbstractVector && length(u) == net.N ||
        throw(ArgumentError("Corrientes inválidas."))
    all(isapprox(u[i], u[1]; atol=1e-12, rtol=0) for i in 2:net.N) ||
        throw(ArgumentError("La cota del cuaderno requiere corriente común."))
    lambda = _law_at(net.lambda, t)
    N = net.N
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    return (1+lambda)*sum((x[2i-1]-vm)^2 for i in 1:N) -
           net.neuron.b*sum((x[2i]-wm)^2 for i in 1:N)
end
