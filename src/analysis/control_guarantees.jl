# ETAPA D — Garantías suficientes y cotas certificadas por hipótesis.

"""    sufficient_decay_rate(p; lambda, k=0)

Si b>0 y epsilon=k-1-lambda>0, devuelve c=min(epsilon,b/tau).
Para control sólo en lambda usar k=0; para feedback usar k>0.
Si las hipótesis no se cumplen se lanza ArgumentError (no inventa garantía).
"""
function sufficient_decay_rate(p::FHNParameters; lambda::Real, k::Real=0.0)
    isfinite(lambda) && isfinite(k) && k >= 0 ||
        throw(ArgumentError("lambda y k deben ser finitos, k>=0."))
    p.b > 0 || throw(ArgumentError("La demostración exige b>0."))
    eps = k - 1 - lambda
    eps > 0 ||
        throw(ArgumentError("No se satisface k > 1+lambda; no hay garantía global."))
    return min(eps, p.b/p.tau)
end

"""    sufficient_arrival_time(E0, delta, p; lambda, k=0)

Tiempo SUFICIENTE del cuaderno:
max{0, log(C_tau*E0/delta)/(2c)}, C_tau=max(1,tau)/min(1,tau).
No es tiempo mínimo ni asegura fases geométricas definidas.
"""
function sufficient_arrival_time(E0::Real, delta::Real, p::FHNParameters;
        lambda::Real, k::Real=0.0)
    isfinite(E0) && E0 >= 0 && isfinite(delta) && delta > 0 ||
        throw(ArgumentError("E0 debe ser >=0 y delta >0, ambos finitos."))
    c = sufficient_decay_rate(p; lambda=lambda, k=k)
    E0 <= delta && return 0.0
    Ctau = max(1.0, p.tau) / min(1.0, p.tau)
    return max(0.0, log(Ctau*E0/delta)/(2c))
end

function _common_baseline(net::FHNNetwork, t::Real)
    u = _current_at(net.currents, t)
    u isa AbstractVector && length(u) == net.N && all(isfinite,u) ||
        return false
    return all(isapprox(ui, u[1]; atol=1e-12, rtol=0) for ui in u)
end

"""    feedback_lyapunov_derivative(x, system; t=0.0)

Diferencia exacta dV/dt usando el campo controlado, sin diferencia finita.
"""
function feedback_lyapunov_derivative(x::AbstractVector,
        system::FHNCurrentFeedback; t::Real=0.0)
    N = system.base.N
    length(x) == 2N || throw(DimensionMismatch("Estado incorrecto."))
    dx = feedback_rhs(x, system, t)
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    tau = system.base.neuron.tau
    return sum((x[2i-1]-vm)*dx[2i-1] +
               tau*(x[2i]-wm)*dx[2i] for i in 1:N)
end

"""    feedback_lyapunov_upper_bound(x, system; t=0.0)

V̇ <= (1+lambda-k) sum_i e_i² - b sum_i h_i².
Sólo se usa bajo campo medio y corriente basal común.
"""
function feedback_lyapunov_upper_bound(x::AbstractVector,
        system::FHNCurrentFeedback; t::Real=0.0)
    net = system.base
    is_mean_field(net) && _common_baseline(net,t) ||
        throw(ArgumentError("La cota exige campo medio y corriente basal común."))
    N = net.N
    length(x) == 2N || throw(DimensionMismatch("Estado incorrecto."))
    vm = sum(x[1:2:end])/N
    wm = sum(x[2:2:end])/N
    lambda = _law_at(net.lambda,t)
    return (1+lambda-system.k)*sum((x[2i-1]-vm)^2 for i in 1:N) -
           net.neuron.b*sum((x[2i]-wm)^2 for i in 1:N)
end

"""    lyapunov_certificate(sol, net_or_feedback)

Devuelve V(t), E(t) y las cotas garantizadas para TODOS los tiempos,
cuando el acoplamiento es constante O de la familia acotada
`AttractiveModulatedCoupling`, la red es todos-con-todos y las corrientes
basales son comunes. No certifica funciones temporales arbitrarias por
comprobar sólo unos pocos puntos.
"""
function lyapunov_certificate(sol, system::Union{FHNNetwork,FHNCurrentFeedback})
    net = system isa FHNNetwork ? system : system.base
    gain = system isa FHNNetwork ? 0.0 : system.k
    is_mean_field(net) ||
        throw(ArgumentError("La garantía exige conectividad todos-con-todos."))
    (net.lambda isa Real || net.lambda isa AttractiveModulatedCoupling) ||
        throw(ArgumentError("El certificado exige lambda constante o modulación con cota global analítica."))
    net.currents isa AbstractVector && _common_baseline(net,first(sol.t)) ||
        throw(ArgumentError("El certificado automático exige corriente basal constante y común; no certifica funciones temporales arbitrarias."))
    p = net.neuron
    upper_lambda = net.lambda isa Real ? net.lambda : net.lambda.ceiling
    c = sufficient_decay_rate(p; lambda=upper_lambda, k=gain)
    V = Float64[lyapunov_dispersion(x,p.tau) for x in sol.u]
    E = Float64[synchronization_error(x) for x in sol.u]
    time = Float64.(sol.t)
    V0, E0 = first(V), first(E)
    Ctau = max(1.0,p.tau)/min(1.0,p.tau)
    decay = exp.(-2c.*(time .- first(time)))
    Vupper = V0 .* decay
    Eupper = (Ctau*E0) .* decay
    return (t=time, c=c, epsilon=gain-1-upper_lambda,
            certified_lambda_upper=upper_lambda,
            V=V, E=E, V_upper=Vupper, E_upper=Eupper,
            Ctau=Ctau,
            max_V_violation=maximum(V .- Vupper),
            max_E_violation=maximum(E .- Eupper))
end
