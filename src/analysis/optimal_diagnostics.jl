# Etapa E - verificación continua por EDO, reconstrucción de adjuntos y shooting.
# Las pruebas de PMP son condiciones NECESARIAS; no certifican óptimo global.

function _ocp_tstops(T,M)
    return [T*j/M for j in 1:M-1]
end

"""Integración independiente con Tsit5 de un control de acoplamiento a trozos."""
function simulate_optimal_sync(s::SyncOCPSpec,result;
        saveat::Real=0.2,reltol::Real=1e-9,abstol::Real=1e-9)
    length(result.controls)==result.M || throw(DimensionMismatch("Control de longitud incorrecta."))
    T=result.T
    law=t -> piecewise_value(result.controls,t,T)
    net=FHNNetwork(s.N,s.neuron;lambda=law)
    problem=ODEProblem(network_rhs!,copy(s.x0),(0.0,T),net)
    sol=solve(problem,Tsit5();saveat=saveat,tstops=_ocp_tstops(T,result.M),
              reltol=reltol,abstol=abstol)
    return (sol=sol,network=net,metrics=network_metrics(sol,net;center=s.center))
end

"""Integración independiente con Tsit5 de corrientes a trozos y suma cero."""
function simulate_optimal_desync(s::DesyncOCPSpec,result;
        saveat::Real=0.2,reltol::Real=1e-9,abstol::Real=1e-9)
    size(result.controls)==(s.N,result.M) || throw(DimensionMismatch("Control N×M incorrecto."))
    law=t -> collect(piecewise_currents(result.controls,t,s.T))
    net=FHNNetwork(s.N,s.neuron;lambda=s.lambda0,currents=law)
    problem=ODEProblem(network_rhs!,copy(s.x0),(0.0,s.T),net)
    sol=solve(problem,Tsit5();saveat=saveat,tstops=_ocp_tstops(s.T,result.M),
              reltol=reltol,abstol=abstol)
    return (sol=sol,network=net,metrics=network_metrics(sol,net;center=s.center))
end

function _solve_adjoint(sol,terminal_adjoint,rhs,T,M;reltol=1e-9,abstol=1e-9)
    prob=ODEProblem((dp,p,par,t)->(dp .= rhs(sol(t),p,t)),
                    Float64.(terminal_adjoint),(T,0.0))
    stops=reverse(_ocp_tstops(T,M))
    return solve(prob,Tsit5();tstops=stops,reltol=reltol,abstol=abstol)
end

"""Diagnóstico PMP para trayectoria de sincronización candidata.

mu proviene de la solución aproximada del problema discretizado.
Se reporta H(T) y residuo de ley proyectada, sin forzar su anulación.
"""
function pmp_diagnostics_sync(s::SyncOCPSpec,result,sol;
        samples::Integer=101,reltol=1e-9,abstol=1e-9)
    samples>=3 || throw(ArgumentError("samples>=3."))
    T=result.T;M=result.M;mu=result.mu
    terminal=mu .* sync_error_gradient(sol.u[end])
    law=t->piecewise_value(result.controls,t,T)
    psol=_solve_adjoint(sol,terminal,
        (x,p,t)->sync_adjoint(x,p,law(t),s),T,M;reltol=reltol,abstol=abstol)
    ts=collect(range(0.0,T,length=samples))
    λ=Float64[law(t) for t in ts]
    λpmp=Float64[sync_control_pmp(sol(t),psol(t),s) for t in ts]
    hh=Float64[sync_hamiltonian(sol(t),psol(t),λ[j],s) for (j,t) in enumerate(ts)]
    pT=psol(T)
    return (kind=:synchronization,t=ts,lambda=λ,lambda_pmp=λpmp,
            H=hh,terminal_E=Float64(synchronization_error(sol.u[end])),
            terminal_error=Float64(synchronization_error(sol.u[end])-s.delta),
            terminal_adjoint_residual=norm(pT.-mu.*sync_error_gradient(sol.u[end])),
            terminal_H=hh[end],max_projected_control_residual=maximum(abs.(λ.-λpmp)),
            mu=mu,adjoint_solution=psol)
end

"""Diagnóstico PMP de desincronización; incluye condición terminal κ∇r²."""
function pmp_diagnostics_desync(s::DesyncOCPSpec,result,sol;
        samples::Integer=101,reltol=1e-9,abstol=1e-9)
    samples>=3 || throw(ArgumentError("samples>=3."))
    pT=s.kappa .* phase_order_squared_gradient(sol.u[end],s.center)
    law=t -> piecewise_currents(result.controls,t,s.T)
    psol=_solve_adjoint(sol,pT,
        (x,p,t)->desync_adjoint(x,p,law(t),s),s.T,result.M;
        reltol=reltol,abstol=abstol)
    ts=collect(range(0.0,s.T,length=samples))
    U=reduce(hcat,[Float64.(law(t)) for t in ts])
    Upmp=reduce(hcat,[desync_control_pmp(psol(t),s) for t in ts])
    H=Float64[desync_hamiltonian(sol(t),psol(t),view(U,:,j),s)
              for (j,t) in enumerate(ts)]
    return (kind=:desynchronization,t=ts,U=U,U_pmp=Upmp,H=H,
        terminal_r2=Float64(_ocp_r2(sol.u[end],s.center)),
        terminal_adjoint_residual=norm(psol(s.T)-pT),
        max_projected_control_residual=maximum(abs.(U.-Upmp)),
        max_current_sum=maximum(abs.(vec(sum(U;dims=1)))),
        max_current=maximum(abs,U),adjoint_solution=psol)
end

# Shooting directo de las ecuaciones necesarias normales en dimensión pequeña.
# No se usa como prueba de optimalidad ni se supone que el Newton converja.
function _sync_pmp_shoot_residual(y,s::SyncOCPSpec)
    n=2*s.N
    T=exp(y[n+1]);mu=exp(y[n+2])
    initial=vcat(s.x0 .+ zero(y[1]), y[1:n])
    function system!(dw,w,_,q)
        x=view(w,1:n);p=view(w,n+1:2n)
        lam=sync_control_pmp(x,p,s)
        dw[1:n] .= T.*sync_ocp_field(x,lam,s.neuron)
        dw[n+1:2n] .= T.*sync_adjoint(x,p,lam,s)
        return nothing
    end
    prob=ODEProblem(system!,initial,(0.0,1.0))
    sol=solve(prob,Tsit5();reltol=1e-8,abstol=1e-8,save_everystep=false)
    xt=view(sol.u[end],1:n);pt=view(sol.u[end],n+1:2n)
    lamT=sync_control_pmp(xt,pt,s)
    return vcat(pt.-mu.*sync_error_gradient(xt),
                _ocp_E(xt)-s.delta,
                sync_hamiltonian(xt,pt,lamT,s))
end

"""    shoot_sync_pmp(spec; p0=zeros(2N), T0=30, mu0=1, ...)

Shooting indirecto opcional para N pequeño, con 2N+2 incógnitas.
La solución existe sólo si el solver converge y los residuos son pequeños.
Para iniciar: p0=diag.adjoint_solution(0), T0=result.T, mu0=result.mu.
"""
function shoot_sync_pmp(s::SyncOCPSpec;
        p0=zeros(2*s.N),T0=30.0,mu0=1.0,maxiters=35,tol=1e-5)
    length(p0)==2*s.N || throw(DimensionMismatch("p0 requiere 2N componentes."))
    all(isfinite,p0) && 0<T0 && 0<mu0 && isfinite(T0) && isfinite(mu0) ||
        throw(ArgumentError("Semillas de shooting inválidas."))
    guess=vcat(Float64.(p0),log(T0),log(mu0))
    problem=NonlinearSolve.NonlinearProblem((y,par)->_sync_pmp_shoot_residual(y,s),guess)
    try
        root=NonlinearSolve.solve(problem,NonlinearSolve.NewtonRaphson();
                                 abstol=tol,reltol=tol,maxiters=maxiters)
        r=_sync_pmp_shoot_residual(root.u,s)
        err=norm(r,Inf)
        return (status=err<=tol ? :converged : :not_converged,
                residual=Float64.(r),residual_norm=Float64(err),
                p0=Float64.(root.u[1:2*s.N]),T=Float64(exp(root.u[end-1])),
                mu=Float64(exp(root.u[end])),retcode=string(root.retcode))
    catch err
        return (status=:solver_error,reason=sprint(showerror,err),
                residual_norm=Inf,p0=Float64.(p0),T=Float64(T0),mu=Float64(mu0))
    end
end
