# Etapa B: cálculo de Hopf, detección de oscilaciones y shooting periódico.
# Depende solamente de los módulos de la etapa A y OrdinaryDiffEq/LinearAlgebra.

"""    hopf_analysis(p; tol=1e-12)

Analiza ambas bifurcaciones de Hopf para `0 < b < 1`, `tau > b`.
Devuelve `nothing` si esas hipótesis no se satisfacen. La clasificación usa el
coeficiente cúbico radial de la forma normal del cuaderno teórico, no un ajuste
basado en simulaciones.
"""
function hopf_analysis(p::FHNParameters; tol::Real=1e-12)
    tol >= 0 && isfinite(tol) || throw(ArgumentError("tol debe ser finita y no negativa"))
    h = hopf_currents(p)
    h === nothing && return nothing
    b, tau = p.b, p.tau
    kappa = ((2b - 1) * tau - b^2) / (8 * (tau - b^2))
    kind = kappa < -tol ? :supercritical : kappa > tol ? :subcritical : :degenerate
    s = sqrt(1 - b/tau)
    slope(v) = -v * b / (1 - b + b*v^2)  # d(Re eigenvalue)/dI en Hopf
    return (I_minus=h.I_minus, I_plus=h.I_plus,
            v_minus=-s, v_plus=s, omega=h.omega,
            period_linear=2pi/h.omega, kappa=kappa, kind=kind,
            slope_minus=slope(-s), slope_plus=slope(s))
end

"""    upward_crossings(sol, level; after=0.0)

Tiempos de cruce de `v(t) = level` en sentido creciente. Se usa interpolación
densa de la solución de SciML para estimar además los estados de cruce.
Los cruces muy tangenciales se ignoran. No supone periodicidad.
"""
function upward_crossings(sol, level::Real; after::Real=0.0)
    times = Float64[]
    states = Vector{Float64}[]
    t = sol.t
    y = sol.u
    for k in 1:(length(t)-1)
        t[k+1] >= after || continue
        v0, v1 = y[k][1] - level, y[k+1][1] - level
        if v0 < 0 && v1 >= 0 && v1 > v0
            a, b = t[k], t[k+1]
            # Bisección con la interpolación continua de la solución.
            for _ in 1:36
                c = (a+b)/2
                if sol(c)[1] < level
                    a = c
                else
                    b = c
                end
            end
            tc = (a+b)/2
            z = Float64.(sol(tc))
            push!(times, tc)
            push!(states, z)
        end
    end
    return (times=times, states=states)
end

"""    detect_limit_cycle(p; u0=(-1.4,-0.2), transient=500.0, window=450.0,
                              dt=0.25, min_crossings=5, period_rtol=0.025,
                              amplitude_rtol=0.06, min_amplitude=0.04)

Integra, descarta un transitorio, toma cruces ascendentes por `v = v*` y
comprueba estabilidad del período y de la amplitud por vueltas sucesivas.
Devuelve un *candidato periódico* si ambas verificaciones pasan, no una prueba
matemática de existencia. Para certificar numéricamente su cierre usar
`refine_periodic_orbit`.
"""
function detect_limit_cycle(p::FHNParameters;
        u0=(-1.4, -0.2), transient::Real=500.0, window::Real=450.0,
        dt::Real=0.25, min_crossings::Integer=5,
        period_rtol::Real=0.025, amplitude_rtol::Real=0.06,
        min_amplitude::Real=0.04,
        reltol::Real=1e-9, abstol::Real=1e-9)
    transient >= 0 && window > 0 && dt > 0 && min_crossings >= 4 ||
        throw(ArgumentError("Se requieren transient>=0, window>0, dt>0 y min_crossings>=4"))
    period_rtol > 0 && amplitude_rtol > 0 && min_amplitude >= 0 ||
        throw(ArgumentError("Tolerancias inválidas"))
    section = equilibrium(p)[1]
    sol = simulate_individual(p; u0=u0, tspan=(0.0, transient+window),
                              saveat=dt, reltol=reltol, abstol=abstol)
    cr = upward_crossings(sol, section; after=transient)
    found = length(cr.times)
    base = (status=:inconclusive, current=p.I, section=section,
            crossings=found, period=NaN, amplitude_v=NaN,
            amplitude_w=NaN, relative_period_spread=NaN,
            relative_amplitude_spread=NaN, initial_state=Float64[],
            crossing_time=NaN, solution=sol)
    if found < min_crossings
        return merge(base, (status=:insufficient_crossings,))
    end

    # Usamos los últimos tres períodos COMPLETOS, no el transitorio.
    ct = cr.times[end-3:end]
    periods = diff(ct)
    period = sum(periods)/length(periods)
    period_spread = maximum(abs.(periods .- period))/period
    amps_v = Float64[]
    amps_w = Float64[]
    for j in 1:3
        # Muestra la órbita con interpolación densa para minimizar sesgo de dt.
        ts = range(ct[j], ct[j+1], length=301)
        z = [sol(t) for t in ts]
        push!(amps_v, maximum(x[1] for x in z) - minimum(x[1] for x in z))
        push!(amps_w, maximum(x[2] for x in z) - minimum(x[2] for x in z))
    end
    av = sum(amps_v)/length(amps_v)
    aw = sum(amps_w)/length(amps_w)
    amplitude_spread = maximum(abs.(amps_v .- av))/max(av, eps())
    # Evita falsos positivos de transitorios o soluciones que decaen al equilibrio.
    accepted = isfinite(period) && av >= min_amplitude &&
               period_spread <= period_rtol &&
               amplitude_spread <= amplitude_rtol
    status = accepted ? :periodic_candidate : :not_converged
    return merge(base,
        (status=status, period=period, amplitude_v=av, amplitude_w=aw,
         relative_period_spread=period_spread,
         relative_amplitude_spread=amplitude_spread,
         initial_state=cr.states[end], crossing_time=cr.times[end]))
end

# Flujo y matriz fundamental: y=(v,w,Phi_11,Phi_21,Phi_12,Phi_22).
# El Jacobiano de flujo en el sistema de shooting se calcula de forma analítica
# mediante la ecuación variacional, no diferenciando el integrador.
function _variational_rhs!(dy, y, p::FHNParameters, t)
    dv, dw = fhn_vector_field(y[1], y[2], p)
    dy[1] = dv
    dy[2] = dw
    a11 = 1 - y[1]^2
    a21 = 1 / p.tau
    a22 = -p.b / p.tau
    dy[3] = a11*y[3] - y[4]
    dy[4] = a21*y[3] + a22*y[4]
    dy[5] = a11*y[5] - y[6]
    dy[6] = a21*y[5] + a22*y[6]
    return nothing
end

function _flow_and_monodromy(p::FHNParameters, z0::AbstractVector, T::Real;
                              reltol=1e-10, abstol=1e-10)
    y0 = [Float64(z0[1]), Float64(z0[2]), 1.0, 0.0, 0.0, 1.0]
    prob = ODEProblem(_variational_rhs!, y0, (0.0, Float64(T)), p)
    sol = solve(prob, Tsit5(); reltol=reltol, abstol=abstol,
                save_everystep=false, save_start=false, save_end=true)
    yf = sol.u[end]
    return (z=yf[1:2], phi=[yf[3] yf[5]; yf[4] yf[6]])
end

function _shooting_residual(p::FHNParameters, y::AbstractVector, section::Real;
                            reltol=1e-10, abstol=1e-10)
    z0 = y[1:2]
    f = _flow_and_monodromy(p, z0, y[3]; reltol=reltol, abstol=abstol)
    return ([f.z[1]-z0[1], f.z[2]-z0[2], z0[1]-section], f)
end

"""    refine_periodic_orbit(p, initial_state, period; section=equilibrium(p)[1], ...)

Newton por *single shooting* con condición de fase `v(0)=section` y sistema
variacional. Resuelve `Phi_T(z)-z=0` y rechaza soluciones de equilibrio.
La convergencia se informa mediante `status`, `residual` y multiplicadores
Floquet; no se presume convergencia si Newton falla.
"""
function refine_periodic_orbit(p::FHNParameters, initial_state, period::Real;
        section::Real=equilibrium(p)[1], tol::Real=1e-8,
        maxiter::Integer=14, min_amplitude::Real=0.025,
        reltol::Real=1e-10, abstol::Real=1e-10)
    length(initial_state) == 2 || throw(ArgumentError("Estado inicial debe ser [v,w]"))
    isfinite(period) && period > 0 && isfinite(section) ||
        throw(ArgumentError("Período positivo y sección finita requeridos"))
    maxiter > 0 && tol > 0 && min_amplitude >= 0 ||
        throw(ArgumentError("Parámetros de Newton inválidos"))
    z0 = Float64.(collect(initial_state))
    all(isfinite, z0) || throw(ArgumentError("Estado inicial no finito"))
    # Ajustar v0 a la sección. El cruce ascendente proviene normalmente de detect_limit_cycle.
    y = [Float64(section), z0[2], Float64(period)]
    status = :max_iterations
    residual = Inf
    iterations = 0
    minT = max(0.1, 0.15*period)
    for k in 1:maxiter
        iterations = k
        r, flow = _shooting_residual(p, y, section; reltol=reltol, abstol=abstol)
        residual = norm(r)
        if residual <= tol
            status = :converged
            break
        end
        J = zeros(3,3)
        J[1:2,1:2] .= flow.phi - Matrix{Float64}(I,2,2)
        J[1:2,3] .= fhn_rhs(flow.z, p)
        J[3,1] = 1.0
        if !all(isfinite, J) || abs(det(J)) < 1e-13
            status = :singular_jacobian
            break
        end
        step = try
            J \ r
        catch
            status = :singular_jacobian
            break
        end
        improved = false
        for i in 0:9
            eta = 2.0^(-i)
            trial = y .- eta .* step
            if trial[3] <= minT || !all(isfinite, trial)
                continue
            end
            newr, _ = _shooting_residual(p, trial, section;
                                         reltol=reltol, abstol=abstol)
            if all(isfinite, newr) && norm(newr) < residual
                y = trial
                improved = true
                break
            end
        end
        if !improved
            status = :no_decrease
            break
        end
    end
    if status != :converged
        r, _ = _shooting_residual(p, y, section; reltol=reltol, abstol=abstol)
        residual = norm(r)
        if residual <= tol
            status = :converged
        end
    end
    # Convergencia a una trayectoria constante no es una órbita periódica.
    zstar = equilibrium(p)
    if status == :converged && norm(y[1:2] .- zstar) < min_amplitude/10
        status = :trivial_equilibrium
    end
    if status != :converged
        return (status=status, iterations=iterations, residual=residual,
                current=p.I, period=y[3], initial_state=copy(y[1:2]),
                amplitude_v=NaN, amplitude_w=NaN,
                floquet=ComplexF64[], phase_multiplier=NaN+0im,
                transverse_multiplier=NaN+0im, stability=:unknown,
                monodromy=zeros(2,2), solution=nothing)
    end
    orbit = simulate_individual(p; u0=y[1:2], tspan=(0.0,y[3]),
                                saveat=y[3]/400, reltol=reltol, abstol=abstol)
    ampv = maximum(u[1] for u in orbit.u) - minimum(u[1] for u in orbit.u)
    ampw = maximum(u[2] for u in orbit.u) - minimum(u[2] for u in orbit.u)
    if ampv < min_amplitude
        return (status=:trivial_equilibrium, iterations=iterations,
                residual=residual, current=p.I, period=y[3],
                initial_state=copy(y[1:2]), amplitude_v=ampv,
                amplitude_w=ampw, floquet=ComplexF64[],
                phase_multiplier=NaN+0im, transverse_multiplier=NaN+0im,
                stability=:unknown, monodromy=zeros(2,2), solution=nothing)
    end
    flow = _flow_and_monodromy(p, y[1:2], y[3]; reltol=reltol, abstol=abstol)
    multipliers = ComplexF64.(eigvals(flow.phi))
    phase_id = argmin(abs.(multipliers .- 1))
    other_id = 3 - phase_id
    phase = multipliers[phase_id]
    transverse = multipliers[other_id]
    stability = abs(phase - 1) > 5e-3 ? :undetermined :
                abs(transverse) < 1 - 1e-4 ? :stable :
                abs(transverse) > 1 + 1e-4 ? :unstable : :near_critical
    return (status=:converged, iterations=iterations, residual=residual,
            current=p.I, period=y[3], initial_state=copy(y[1:2]),
            amplitude_v=ampv, amplitude_w=ampw, floquet=multipliers,
            phase_multiplier=phase, transverse_multiplier=transverse,
            stability=stability, monodromy=flow.phi, solution=orbit)
end

"""    refine_periodic_orbit(p, detection::NamedTuple; kwargs...)

Refina un resultado de `detect_limit_cycle` solo si los chequeos pasaron.
"""
function refine_periodic_orbit(p::FHNParameters, detection::NamedTuple; kwargs...)
    detection.status === :periodic_candidate ||
        throw(ArgumentError("No hay candidato periódico confirmado por la detección"))
    return refine_periodic_orbit(p, detection.initial_state, detection.period;
                                  section=detection.section, kwargs...)
end

"""    continue_periodic_branch(p, orbit, currents; tol=1e-8)

Continuación por parámetros (pasos en I) de una órbita obtenida por shooting.
Usa el último estado periódico como predictor y vuelve a corregir por Newton.
No es continuación de pseudo-arclength: cerca de pliegues de ciclos puede fallar.
Cada elemento reporta su estado y precisión; la lista se detiene en el primer fallo.
"""
function continue_periodic_branch(p::FHNParameters, orbit, currents;
                                  tol::Real=1e-8, maxiter::Integer=15)
    orbit.status == :converged ||
        throw(ArgumentError("Se requiere una órbita inicial refinada"))
    records = NamedTuple[]
    current_orbit = orbit
    for Ival in currents
        q = with_current(p, Ival)
        predicted = copy(current_orbit.initial_state)
        predicted[1] = equilibrium(q)[1]
        nxt = refine_periodic_orbit(q, predicted, current_orbit.period;
                                    tol=tol, maxiter=maxiter)
        push!(records, (I=q.I, orbit=nxt))
        nxt.status === :converged || break
        current_orbit = nxt
    end
    return records
end

"""    hopf_seed(p, side; offset=0.003)

Predicción local a partir de la forma normal radial cerca de uno de los Hopf.
Útil para iniciar shooting de órbitas pequeñas (inclusive inestables).
NO asegura convergencia del método indirecto ni reemplaza continuación robusta.
"""
function hopf_seed(p::FHNParameters, side::Symbol; offset::Real=0.003)
    h = hopf_analysis(p)
    h === nothing && throw(ArgumentError("No hay candidatos de Hopf no degenerados"))
    h.kind == :degenerate && throw(ArgumentError("Hopf degenerado: falta orden superior"))
    side in (:minus, :plus) || throw(ArgumentError("side debe ser :minus o :plus"))
    isfinite(offset) && offset > 0 || throw(ArgumentError("offset debe ser positivo"))
    Ih = side == :minus ? h.I_minus : h.I_plus
    slope = side == :minus ? h.slope_minus : h.slope_plus
    # La forma normal requiere sigma * kappa < 0.
    direction = sign(-h.kappa / slope)
    q = with_current(p, Ih + direction * offset)
    sigma = stability_info(q).trace / 2
    sigma * h.kappa < 0 || throw(ArgumentError("Punto fuera del lado con ciclo pequeño"))
    radius = sqrt(-sigma / h.kappa)
    zstar = equilibrium(q)
    # En X=0,Y=-R: desplazamiento w=-omega*R; v crece al cruzar sección.
    return (params=q, side=side, initial_state=[zstar[1], zstar[2]-h.omega*radius],
            period=h.period_linear, amplitude_estimate=2radius)
end
