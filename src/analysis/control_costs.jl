# ETAPA D — Evaluación (NO optimización) de las funciones objetivo del PDF.

function _trapezoidal_integral(t::AbstractVector, f::AbstractVector)
    length(t) == length(f) && length(t) >= 1 ||
        throw(DimensionMismatch("Tiempo e integrando incompatibles."))
    all(isfinite, t) && all(isfinite, f) ||
        throw(ArgumentError("Datos no finitos en la cuadratura."))
    total = 0.0
    for j in 2:length(t)
        dt = t[j] - t[j-1]
        dt > 0 || throw(ArgumentError("Los tiempos deben ser crecientes."))
        total += dt*(f[j-1]+f[j])/2
    end
    return total
end

function _cumulative_trapezoid(t::AbstractVector, f::AbstractVector)
    length(t) == length(f) && length(t) >= 1 ||
        throw(DimensionMismatch("Tiempo e integrando incompatibles."))
    acc = zeros(Float64, length(t))
    for j in 2:length(t)
        dt = t[j]-t[j-1]
        dt > 0 || throw(ArgumentError("Los tiempos deben ser crecientes."))
        acc[j] = acc[j-1] + dt*(f[j-1]+f[j])/2
    end
    return acc
end

"""    first_tolerance_hit(sol, delta; max_iter=50)

Primera detección de E(x(t)) <= delta usando puntos guardados de la solución;
refina el cruce mediante bisección y la solución continua `sol(t)`.
Si ocurre un cruce completamente entre dos muestras y vuelve a salir, puede
no detectarse: elegir saveat suficientemente pequeño.
"""
function first_tolerance_hit(sol, delta::Real; max_iter::Integer=50)
    isfinite(delta) && delta > 0 ||
        throw(ArgumentError("delta debe ser finita y >0."))
    max_iter > 0 || throw(ArgumentError("max_iter debe ser positivo."))
    n = length(sol.t)
    n > 0 || throw(ArgumentError("Solución vacía."))
    if synchronization_error(sol.u[1]) <= delta
        return (status=:already_reached, time=Float64(sol.t[1]),
                state=Float64.(sol.u[1]), E=synchronization_error(sol.u[1]),
                saved_index=1)
    end
    for j in 2:n
        if synchronization_error(sol.u[j]) <= delta
            a, b = Float64(sol.t[j-1]), Float64(sol.t[j])
            for _ in 1:max_iter
                mid = (a+b)/2
                if synchronization_error(sol(mid)) > delta
                    a = mid
                else
                    b = mid
                end
            end
            hit_time = (a+b)/2
            hit_state = Float64.(sol(hit_time))
            return (status=:reached, time=hit_time,
                    state=hit_state, E=synchronization_error(hit_state),
                    saved_index=j)
        end
    end
    return (status=:not_reached, time=NaN, state=nothing,
            E=Float64(synchronization_error(sol.u[end])), saved_index=0)
end

function _check_positive_cost_weights(values)
    all(x -> x isa Real && isfinite(x) && x > 0, values) ||
        throw(ArgumentError("Todos los pesos deben ser reales, finitos y positivos."))
end

function _sync_running_parts(x, t, net::FHNNetwork, center,
        rho, alpha, beta, gamma, lambda_ref)
    r2 = phase_order_squared(x, center)
    isfinite(r2) ||
        throw(DomainError(t, "No hay fase definida en este instante; J no es evaluable."))
    lambda = _law_at(net.lambda,t)
    isfinite(lambda) || throw(ArgumentError("lambda(t) debe ser finita."))
    return (time=Float64(rho),
            state=Float64(alpha*synchronization_error(x)),
            phase=Float64(beta*(1-r2)^2),
            effort=Float64(gamma*(lambda-lambda_ref)^2))
end

"""    evaluate_sync_cost(sol, net; delta, rho=1, alpha=0.1, beta=0.2,
                              gamma=0.5, lambda_ref=-0.15, center=...)

Evalúa J_s = ∫[rho+alpha E+beta(1-r²)²+gamma(lambda-lambda_ref)²]dt.
- `horizon_J` evalúa hasta el fin de la simulación (diagnóstico comparativo).
- `arrival_J` integra sólo hasta el PRIMER E=delta (T libre), o Inf si no llega.
- No minimiza J ni afirma optimalidad. Exige corrientes COMUNES.
La integral numérica se calcula por trapecios, incorporando el cruce refinado.
"""
function evaluate_sync_cost(sol, net::FHNNetwork;
        delta::Real=0.01, rho::Real=1.0, alpha::Real=0.1,
        beta::Real=0.2, gamma::Real=0.5, lambda_ref::Real=-0.15,
        center=equilibrium(net.neuron))
    _check_positive_cost_weights((rho,alpha,beta,gamma))
    isfinite(lambda_ref) || throw(ArgumentError("lambda_ref debe ser finita."))
    isfinite(delta) && delta > 0 || throw(ArgumentError("delta debe ser >0."))
    all(t -> _common_baseline(net,t), sol.t) ||
        throw(ArgumentError("Psin requiere corriente común, no corrientes diferenciadas."))
    length(center) == 2 && all(isfinite,center) ||
        throw(ArgumentError("Centro de fase inválido."))
    times = Float64.(sol.t)
    parts = [_sync_running_parts(sol.u[j],times[j],net,center,
             rho,alpha,beta,gamma,lambda_ref) for j in eachindex(times)]
    L = Float64[sum(values(p)) for p in parts]
    accumulated = _cumulative_trapezoid(times,L)
    hit = first_tolerance_hit(sol,delta)
    if hit.status == :not_reached
        arrival = Inf
        T = NaN
    elseif hit.status == :already_reached
        arrival = 0.0
        T = times[1]
    else
        # Suma de trapecios completos hasta la muestra anterior al cruce.
        j = hit.saved_index
        running = accumulated[j-1]
        end_parts = _sync_running_parts(hit.state,hit.time,net,center,
                    rho,alpha,beta,gamma,lambda_ref)
        Lend = sum(values(end_parts))
        arrival = running + (hit.time-times[j-1])*(L[j-1]+Lend)/2
        T = hit.time
    end
    components = (
        time=_trapezoidal_integral(times,[p.time for p in parts]),
        state=_trapezoidal_integral(times,[p.state for p in parts]),
        phase=_trapezoidal_integral(times,[p.phase for p in parts]),
        effort=_trapezoidal_integral(times,[p.effort for p in parts]))
    return (status=hit.status, reached=hit.status != :not_reached,
            arrival_time=T, arrival_J=Float64(arrival),
            horizon_J=Float64(last(accumulated)),
            t=times, running_cost=L, accumulated_cost=accumulated,
            horizon_components=components,
            delta=Float64(delta), terminal_E=Float64(synchronization_error(sol.u[end])),
            phase_defined=true)
end

"""    evaluate_desync_cost(sol, net; alpha_d=1, gamma_d=0.2,
                                 kappa=2, Umax=1, center=...)

Evalúa J_d = ∫ [alpha_d*r² + gamma_d*Σu_i²]dt + kappa*r²(T),
con horizonte fijo, lambda0 constante, |u_i|<=Umax, Σu_i=0.
Si falla una restricción se lanza error; no lo presenta como admisible.
"""
function evaluate_desync_cost(sol, net::FHNNetwork;
        alpha_d::Real=1.0, gamma_d::Real=0.2,
        kappa::Real=2.0, Umax::Real=1.0,
        center=equilibrium(net.neuron))
    _check_positive_cost_weights((alpha_d,gamma_d,kappa,Umax))
    net.lambda isa Real ||
        throw(ArgumentError("Pdes exige un acoplamiento lambda0 constante."))
    times = Float64.(sol.t)
    U = control_current_matrix(sol,net)
    r2 = Float64[phase_order_squared(x,center) for x in sol.u]
    all(isfinite,r2) ||
        throw(DomainError(NaN,"Fase indefinida durante la trayectoria."))
    for j in eachindex(times)
        abs(sum(view(U,:,j))) <= 1e-8 ||
            throw(ArgumentError("Se viola sum_i u_i=0."))
        maximum(abs,view(U,:,j)) <= Umax+1e-8 ||
            throw(ArgumentError("Se viola la cota |u_i|<=Umax."))
    end
    energy = vec(sum(abs2,U; dims=1))
    L = alpha_d.*r2 .+ gamma_d.*energy
    acc = _cumulative_trapezoid(times,L)
    terminal_cost = kappa*last(r2)
    return (total_J=Float64(last(acc)+terminal_cost),
            running_J=Float64(last(acc)), terminal_J=Float64(terminal_cost),
            phase_running_J=Float64(alpha_d*_trapezoidal_integral(times,r2)),
            current_running_J=Float64(gamma_d*_trapezoidal_integral(times,energy)),
            t=times, running_cost=L, accumulated_cost=acc,
            r2=r2, current_energy=energy, terminal_r2=last(r2),
            terminal_E=synchronization_error(sol.u[end]),
            initial_E=synchronization_error(sol.u[1]),
            max_current=maximum(abs,U), max_current_sum=maximum(abs.(vec(sum(U;dims=1)))))
end
