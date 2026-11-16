# Etapa E - optimización directa AUTOCONTENIDA por control constante a trozos.
# Integración RK4 + diferenciación automática de J + descenso proyectado.
# La solución discretizada es un CANDIDATO; no es una prueba de óptimo global.

"""Control constante a trozos. En los nodos de salto se toma el intervalo derecho."""
function piecewise_value(values::AbstractVector,t::Real,T::Real)
    M=length(values)
    M>0 && T>0 || throw(ArgumentError("Controles no vacíos y T>0 requeridos."))
    idx=clamp(floor(Int,Float64(t)/Float64(T)*M)+1,1,M)
    return values[idx]
end

function piecewise_currents(U::AbstractMatrix,t::Real,T::Real)
    idx=clamp(floor(Int,Float64(t)/Float64(T)*size(U,2))+1,1,size(U,2))
    return view(U,:,idx)
end

# Cuadratura RK4 aplicada también al estado acumulador J, para consistencia.
function _rk4_sync(s::SyncOCPSpec, lambdas::AbstractVector,T,substeps::Int)
    M=length(lambdas)
    h=T/(M*substeps)
    x=s.x0 .+ zero(T+lambdas[1])
    acc=zero(x[1])
    for j in 1:M
        lam=lambdas[j]
        for _ in 1:substeps
            k1=sync_ocp_field(x,lam,s.neuron); c1=sync_running_cost(x,lam,s)
            x2=x .+ (h/2).*k1
            k2=sync_ocp_field(x2,lam,s.neuron); c2=sync_running_cost(x2,lam,s)
            x3=x .+ (h/2).*k2
            k3=sync_ocp_field(x3,lam,s.neuron); c3=sync_running_cost(x3,lam,s)
            x4=x .+ h.*k3
            k4=sync_ocp_field(x4,lam,s.neuron); c4=sync_running_cost(x4,lam,s)
            x=x .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
            acc += h*(c1+2*c2+2*c3+c4)/6
        end
    end
    return (J=acc,terminal_E=_ocp_E(x),terminal_r2=_ocp_r2(x,s.center),xT=x)
end

function _rk4_desync(s::DesyncOCPSpec,U::AbstractMatrix,substeps::Int)
    N,M=size(U)
    N==s.N || throw(DimensionMismatch("U requiere N filas."))
    h=s.T/(M*substeps)
    x=s.x0 .+ zero(U[1,1])
    acc=zero(x[1])
    for j in 1:M
        u=view(U,:,j)
        for _ in 1:substeps
            k1=desync_ocp_field(x,u,s.lambda0,s.neuron);c1=desync_running_cost(x,u,s)
            x2=x .+ (h/2).*k1
            k2=desync_ocp_field(x2,u,s.lambda0,s.neuron);c2=desync_running_cost(x2,u,s)
            x3=x .+ (h/2).*k2
            k3=desync_ocp_field(x3,u,s.lambda0,s.neuron);c3=desync_running_cost(x3,u,s)
            x4=x .+ h.*k3
            k4=desync_ocp_field(x4,u,s.lambda0,s.neuron);c4=desync_running_cost(x4,u,s)
            x=x .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
            acc += h*(c1+2*c2+2*c3+c4)/6
        end
    end
    r2=_ocp_r2(x,s.center)
    return (J=acc+s.kappa*r2, running_J=acc, terminal_J=s.kappa*r2,
            terminal_r2=r2,terminal_E=_ocp_E(x),xT=x)
end

function _check_discretization(M,substeps)
    M>=2 && substeps>=1 || throw(ArgumentError("Se requiere M>=2 y substeps>=1."))
end

function _sync_projection!(y,s::SyncOCPSpec,Tbounds)
    for j in 1:length(y)-1
        y[j]=clamp(y[j],s.lambda_min,s.lambda_max)
    end
    y[end]=clamp(y[end],Tbounds[1],Tbounds[2])
    return y
end

function _desync_projection!(z,N,Umax)
    M=div(length(z),N)
    length(z)==N*M || throw(DimensionMismatch("z debe tener N*M componentes."))
    for j in 1:M
        ids=(j-1)*N+1:j*N
        z[ids] .= project_balanced_box(view(z,ids),Umax)
    end
    return z
end

# Descenso proyectado Armijo sin dependencias adicionales. No declara convergencia
# salvo que el gradiente proyectado sea pequeño; registra todo el historial.
function _projected_descent(f,initial,project!;
        maxiter=70,step0=1.0,grad_tol=5e-5)
    y=copy(Float64.(initial));project!(y)
    fy=f(y)
    best=copy(y);bestval=fy
    history=Float64[fy]
    status=:max_iterations
    for _ in 1:maxiter
        g=ForwardDiff.gradient(f,y)
        all(isfinite,g) || (status=:invalid_gradient;break)
        probe=project!(y .- min(step0,1.0).*g)
        if maximum(abs.(probe.-y))<grad_tol
            status=:projected_stationary
            break
        end
        accepted=false
        step=step0
        for _ in 1:30
            trial=project!(y .- step.*g)
            if maximum(abs.(trial.-y))<1e-14
                break
            end
            ft=try
                f(trial)
            catch err
                if err isa DomainError || err isa OverflowError
                    Inf
                else
                    rethrow(err)
                end
            end
            if isfinite(ft) && ft <= fy + 1e-4*dot(g,trial.-y)
                y=trial;fy=ft;accepted=true
                break
            end
            step *= 0.5
        end
        if !accepted
            status=:line_search_stalled
            break
        end
        if fy < bestval
            bestval=fy;best=copy(y)
        end
        push!(history,fy)
    end
    return (parameters=best,value=bestval,history=history,status=status)
end

# La restricción E(x(T)) <= delta es DURA: un Lagrangiano aumentado
# solo produce CANDIDATOS, nunca convierte una trayectoria infactible en factible.
# El margen evita aceptar como factible una órbita cerca del umbral que luego
# lo supera al integrarla de manera independiente con Tsit5.
_sync_margin(s::SyncOCPSpec) = min(1e-4,0.005*s.delta)

function _sync_eval_candidate(s,z,M,substeps)
    try
        a=_rk4_sync(s,view(z,1:M),z[end],substeps)
        return isfinite(a.J) && isfinite(a.terminal_E) ? a : nothing
    catch err
        if err isa DomainError || err isa OverflowError
            return nothing
        end
        rethrow(err)
    end
end

# Construye controles TESTIGO factibles con lambda constante dentro de las
# cotas. Se busca el primer tiempo muestreado que alcanza la tolerancia,
# y se refina el instante por bisección. No se presupone que E(t) sea monótona.
function _sync_constant_seeds(s::SyncOCPSpec,M,substeps,Tbounds,initial_controls;
                              margin=_sync_margin(s),grid_points=65)
    lambdas=unique([clamp(v,s.lambda_min,s.lambda_max)
                    for v in (sum(initial_controls)/length(initial_controls),
                              -1.5,-1.35,s.lambda_min)])
    seeds=NamedTuple[]
    target=s.delta-margin
    for lam in lambdas
        controls=fill(Float64(lam),M)
        prev_t=nothing
        prev_E=Inf
        for T in range(Tbounds[1],Tbounds[2],length=grid_points)
            y=vcat(controls,Float64(T))
            r=_sync_eval_candidate(s,y,M,substeps)
            if r !== nothing && r.terminal_E<=target
                hi=Float64(T)
                if prev_t !== nothing && prev_E>target
                    lo=Float64(prev_t)
                    for _ in 1:26
                        mid=(lo+hi)/2
                        ym=vcat(controls,mid)
                        rm=_sync_eval_candidate(s,ym,M,substeps)
                        if rm!==nothing && rm.terminal_E<=target
                            hi=mid
                        else
                            lo=mid
                        end
                    end
                end
                yy=vcat(controls,hi)
                rr=_sync_eval_candidate(s,yy,M,substeps)
                if rr !== nothing && rr.terminal_E<=target
                    push!(seeds,(y=yy,J=Float64(rr.J),source=:constant_witness))
                end
                break
            end
            if r !== nothing
                prev_t=Float64(T)
                prev_E=Float64(r.terminal_E)
            end
        end
    end
    return seeds
end

# Recupera factibilidad desde el resultado del Lagrangiano aumentado,
# interpolándolo con controles testigo. No presupone convexidad del conjunto
# factible: se comprueba E(T) para cada punto propuesto.
function _sync_recover_feasibility(s,yaug,witnesses,M,substeps;
                                    margin=_sync_margin(s))
    target=s.delta-margin
    out=NamedTuple[]
    function add_if_feasible!(z,source)
        r=_sync_eval_candidate(s,z,M,substeps)
        if r !== nothing && r.terminal_E<=target
            push!(out,(y=copy(z),J=Float64(r.J),source=source))
            return true
        end
        return false
    end
    add_if_feasible!(yaug,:augmented_lagrangian)
    for seed in witnesses
        add_if_feasible!(seed.y,:constant_witness)
        prev_y=copy(yaug)
        prev_ok=false
        for θ in range(0.0,1.0,length=17)[2:end]
            z=(1-θ).*yaug .+ θ.*seed.y
            now_ok=add_if_feasible!(z,:interpolated_restoration)
            if now_ok && !prev_ok
                # Apretar la restricción por bisección hacia el punto infactible.
                lo=copy(prev_y);hi=copy(z)
                for _ in 1:22
                    mid=(lo.+hi)./2
                    rm=_sync_eval_candidate(s,mid,M,substeps)
                    if rm !== nothing && rm.terminal_E<=target
                        hi=mid
                    else
                        lo=mid
                    end
                end
                add_if_feasible!(hi,:boundary_restoration)
            end
            prev_y=z
            prev_ok=now_ok
        end
    end
    return out
end

# Pequeño refinamiento factible del costo. Rechaza cualquier paso que
# cruce E(T)>delta-margin; no usa penalizaciones como sustituto de restricciones.
function _sync_feasible_refinement(s,z0,M,substeps,Tbounds;
                                   margin=_sync_margin(s),iterations=12)
    z=copy(z0)
    project!(q)=_sync_projection!(q,s,Tbounds)
    rz=_sync_eval_candidate(s,z,M,substeps)
    rz === nothing && return nothing
    rz.terminal_E<=s.delta-margin || return nothing
    Jz=Float64(rz.J)
    for _ in 1:iterations
        gj=ForwardDiff.gradient(q->_rk4_sync(s,view(q,1:M),q[end],substeps).J,z)
        ge=ForwardDiff.gradient(q->_rk4_sync(s,view(q,1:M),q[end],substeps).terminal_E,z)
        all(isfinite,gj) && all(isfinite,ge) || break
        d=.-gj
        if rz.terminal_E>=s.delta-3margin && dot(ge,d)>0
            d .-= (dot(ge,d)/max(dot(ge,ge),eps())) .* ge
        end
        mag=maximum(abs,d)
        mag>1e-12 || break
        d ./= max(1.0,mag)
        accepted=false
        step=1.0
        for _ in 1:28
            trial=project!(z .+ step.*d)
            r=_sync_eval_candidate(s,trial,M,substeps)
            if r !== nothing && r.terminal_E<=s.delta-margin &&
               isfinite(r.J) && r.J < Jz-1e-8
                z=trial;rz=r;Jz=Float64(r.J);accepted=true
                break
            end
            step/=2
        end
        accepted || break
    end
    return (y=z,J=Jz,source=:feasible_refinement)
end

# Verificación INDEPENDIENTE: no basta la misma integración RK4 usada al optimizar.
function _sync_terminal_tsit5(s,z,M)
    r=(controls=copy(z[1:M]),M=M,T=Float64(z[end]))
    sim=simulate_optimal_sync(s,r;saveat=max(0.05,r.T/200),
                              reltol=1e-9,abstol=1e-9)
    return Float64(_ocp_E(sim.sol.u[end]))
end

"""    optimize_sync_direct(spec; ...)

Transcripción directa con tiempo final libre y restricción terminal DURA
`E(x(T)) <= delta`. El Lagrangiano aumentado encuentra un candidato inicial;
si éste es infactible, se recupera factibilidad mediante controles constantes
admisibles e interpolación y se refina sin admitir pasos infactibles.
El resultado sólo declara `feasible=true` cuando RK4 Y Tsit5 satisfacen la
restricción. Si no se encuentra un testigo factible, devuelve `feasible=false`
(con `solver_status=:no_verified_feasible_candidate`). No prueba optimalidad.
"""
function optimize_sync_direct(s::SyncOCPSpec;
        M::Integer=16,substeps::Integer=10,T0::Real=35.0,
        T_bounds=(1.0,120.0),initial_lambda=-1.35,
        penalty0::Real=10.0,outer_iterations::Integer=5,
        inner_iterations::Integer=60,feas_tol::Real=5e-3)
    _check_discretization(M,substeps)
    length(T_bounds)==2 || throw(ArgumentError("T_bounds=(mín,máx)."))
    Tb=Float64.(T_bounds)
    all(isfinite,Tb) && 0<Tb[1]<Tb[2] || throw(ArgumentError("T_bounds incorrectas."))
    isfinite(T0) && Tb[1]<=T0<=Tb[2] || throw(ArgumentError("T0 fuera de cotas."))
    penalty0>0 && outer_iterations>=1 && inner_iterations>=1 && feas_tol>0 ||
        throw(ArgumentError("Parámetros del optimizador incorrectos."))
    controls=initial_lambda isa Real ? fill(Float64(initial_lambda),M) : Float64.(collect(initial_lambda))
    length(controls)==M || throw(DimensionMismatch("initial_lambda requiere M valores."))
    y=vcat(controls,Float64(T0))
    project!(z)=_sync_projection!(z,s,Tb)
    project!(y)
    baseline=_rk4_sync(s,view(y,1:M),y[end],substeps)
    mu=0.0;penalty=Float64(penalty0)
    histories=Float64[];statuses=Symbol[]
    margin=_sync_margin(s)
    for _ in 1:outer_iterations
        mufixed=mu;σfixed=penalty
        function augmented(z)
            zsol=_rk4_sync(s,view(z,1:M),z[end],substeps)
            g=zsol.terminal_E-s.delta
            return zsol.J+(max(zero(g),mufixed+σfixed*g)^2-mufixed^2)/(2σfixed)
        end
        descent=_projected_descent(augmented,y,project!;maxiter=inner_iterations)
        y=descent.parameters
        append!(histories,descent.history)
        push!(statuses,descent.status)
        trial=_rk4_sync(s,view(y,1:M),y[end],substeps)
        g=Float64(trial.terminal_E-s.delta)
        mu=max(0.0,mu+penalty*g)
        if -feas_tol<=g<=-margin && mu>0
            break
        end
        penalty *= 4
    end
    # Nunca escoger por costo un candidato infactible.
    witnesses=_sync_constant_seeds(s,M,substeps,Tb,controls;margin=margin)
    candidates=_sync_recover_feasibility(s,y,witnesses,M,substeps;margin=margin)
    if !isempty(candidates)
        idx=argmin(getproperty.(candidates,:J))
        refined=_sync_feasible_refinement(s,candidates[idx].y,M,substeps,Tb;margin=margin)
        refined!==nothing && push!(candidates,refined)
    end
    sort!(candidates;by=c->c.J)
    chosen=nothing
    independent_E=NaN
    for cand in candidates
        independent_E=try
            _sync_terminal_tsit5(s,cand.y,M)
        catch err
            if err isa DomainError || err isa OverflowError
                NaN
            else
                rethrow(err)
            end
        end
        if isfinite(independent_E) && independent_E<=s.delta
            chosen=cand
            break
        end
    end
    if chosen === nothing
        # Fallo explícito: no se imprime ni guarda como óptimo un resultado inválido.
        final=_rk4_sync(s,view(y,1:M),y[end],substeps)
        source=:no_verified_feasible_candidate
        status=:no_verified_feasible_candidate
        zfinal=copy(y)
    else
        zfinal=chosen.y
        final=_rk4_sync(s,view(zfinal,1:M),zfinal[end],substeps)
        source=chosen.source
        status=source==:augmented_lagrangian ? :feasible_augmented_lagrangian : :feasible_restored
    end
    g=Float64(final.terminal_E-s.delta)
    feasible=(g<=0.0 && isfinite(independent_E) && independent_E<=s.delta && chosen!==nothing)
    # El multiplicador del AL es aproximado; no se ha resuelto necesariamente PMP.
    return (kind=:synchronization,controls=copy(zfinal[1:M]),T=Float64(zfinal[end]),
            M=Int(M),substeps=Int(substeps),J=Float64(final.J),
            initial_J=Float64(baseline.J),terminal_E=Float64(final.terminal_E),
            terminal_E_independent=Float64(independent_E),
            terminal_r2=Float64(final.terminal_r2),state_final=Float64.(final.xT),
            constraint_residual=g,feasible=feasible,mu=mu,
            solver_status=status,feasible_source=source,
            multiplier_status=:approximate_unverified,
            feasibility_margin=margin,
            inner_statuses=statuses,objective_history=histories,
            method=:augmented_lagrangian_with_verified_feasibility_restoration,
            gradient_method=:ForwardDiff,lambda_bounds=(s.lambda_min,s.lambda_max))
end

function _desync_initial_guesses(s::DesyncOCPSpec,M,amplitude)
    N=s.N
    v=zeros(Float64,N,M)
    seeds=[copy(v)]
    for pattern in (N==2 ? (:split,:alternating) : (:split,:alternating,:sinusoidal))
        q=balanced_current_profile(N;pattern=pattern)
        for mode in (1,2)
            B=zeros(Float64,N,M)
            for j in 1:M
                fac=mode==1 ? 1.0 : sin(2pi*(j-0.5)/M)
                B[:,j] .= amplitude*fac.*q
            end
            push!(seeds,B)
        end
    end
    return seeds
end

"""    optimize_desync_direct(spec; M=18, substeps=10, ...)

Minimiza el costo discreto de desincronización con T fijo y las restricciones
u_i en [-Umax,Umax], sum_i u_i=0 EN CADA SEGMENTO. Multistart no nulo evita
quedar atrapado en el extremal trivial p=u=0. Compara también con u=0.
No se garantiza que aparezca desincronización sostenida de fase.
"""
function optimize_desync_direct(s::DesyncOCPSpec;
        M::Integer=18,substeps::Integer=10,
        maxiter::Integer=65,seed_amplitude::Real=0.65,
        additional_seeds=Matrix{Float64}[])
    _check_discretization(M,substeps)
    maxiter>=0 && 0<=seed_amplitude<=s.Umax ||
        throw(ArgumentError("Iteraciones o amplitud inicial inválidas."))
    seeds=_desync_initial_guesses(s,M,seed_amplitude)
    for S in additional_seeds
        size(S)==(s.N,M) || throw(DimensionMismatch("Cada semilla debe ser N×M."))
        push!(seeds,Matrix{Float64}(S))
    end
    project!(z)=_desync_projection!(z,s.N,s.Umax)
    f(z)=_rk4_desync(s,reshape(z,s.N,M),substeps).J
    baseline=f(zeros(s.N*M))
    outcomes=NamedTuple[]
    for (i,seed) in enumerate(seeds)
        result=_projected_descent(f,vec(seed),project!;maxiter=maxiter)
        push!(outcomes,(seed=i,objective=result.value,status=result.status,
                        parameters=result.parameters,history=result.history))
    end
    idx=argmin(getproperty.(outcomes,:objective))
    best=outcomes[idx]
    U=reshape(copy(best.parameters),s.N,M)
    final=_rk4_desync(s,U,substeps)
    maxsum=maximum(abs(sum(U[:,j])) for j in 1:M)
    maxabs=maximum(abs,U)
    return (kind=:desynchronization,controls=U,T=s.T,M=Int(M),
            substeps=Int(substeps),J=Float64(final.J),
            baseline_J=Float64(baseline),running_J=Float64(final.running_J),
            terminal_J=Float64(final.terminal_J),terminal_r2=Float64(final.terminal_r2),
            terminal_E=Float64(final.terminal_E),state_final=Float64.(final.xT),
            max_current=maxabs,max_current_sum=maxsum,
            best_seed=idx,seed_costs=getproperty.(outcomes,:objective),
            solver_status=best.status,objective_history=best.history,
            method=:multistart_projected_RK4,gradient_method=:ForwardDiff)
end
