# HJB dependiente del tiempo para N=2, en su estado COMPLETO de 4 dimensiones.
# Controles admisibles u=(q,-q), |q|<=Umax; paso semi-lagrangiano,
# horizonte corto y dominio compacto. Las fronteras fuera de la grilla se
# consideran NO resueltas (Inf), sin extrapolación arbitraria.

struct SmallDesyncHJB
    spec::DesyncOCPSpec
    axes::NTuple{4,Vector{Float64}}
    actions::Vector{Float64}
    h::Float64
    values::Vector{Array{Float64,4}} # 1 = tiempo 0, end = tiempo T
end

function _small_step(x,q,s::DesyncOCPSpec,h)
    u=[q,-q]
    k1=desync_ocp_field(x,u,s.lambda0,s.neuron)
    xm=x .+ (h/2).*k1
    k2=desync_ocp_field(xm,u,s.lambda0,s.neuron)
    return x .+ h.*k2
end

function _sl_interp4(axes::NTuple{4,Vector{Float64}},W::AbstractArray{<:Real,4},x)
    length(x)==4 || throw(DimensionMismatch("Se requieren 4 coordenadas."))
    k=Vector{Int}(undef,4);t=Vector{Float64}(undef,4)
    for j in 1:4
        a=axes[j]
        (a[1]<=x[j]<=a[end]) || return Inf
        k[j]=clamp(searchsortedlast(a,x[j]),1,length(a)-1)
        t[j]=clamp((x[j]-a[k[j]])/(a[k[j]+1]-a[k[j]]),0.0,1.0)
    end
    acc=0.0
    for corner in 0:15
        indices=ntuple(j -> k[j]+((corner>>(j-1))&1),4)
        weight=1.0
        for j in 1:4
            weight*=((corner>>(j-1))&1)==1 ? t[j] : 1-t[j]
        end
        weight<=1e-14 && continue
        val=W[indices...]
        isfinite(val) || return Inf
        acc+=weight*val
    end
    return acc
end

function _small_cost(x,q,s,h)
    return h*(s.alpha_d*_ocp_r2(x,s.center)+2*s.gamma_d*q^2)
end

function _small_bellman_action(x,q,s,h,axes,nextvalue)
    xn=_small_step(x,q,s,h)
    future=_sl_interp4(axes,nextvalue,xn)
    return isfinite(future) ? _small_cost(x,q,s,h)+future : Inf
end

"""    solve_small_desync_hjb(spec; vaxis, waxis, actions, steps=8)

HJB temporal discreta en la red COMPLETA N=2 (4 variables). La grilla
cartesiana es finita y las corrientes pertenecen a una malla en [-U,U].
No constituye una solución global de viscosidad certificada.
"""
function solve_small_desync_hjb(s::DesyncOCPSpec;
        vaxis=collect(range(-2.4,2.4,length=9)),
        waxis=collect(range(-0.1,2.3,length=7)),
        actions=collect(range(-s.Umax,s.Umax,length=5)),steps=7)
    s.N==2 || throw(ArgumentError("HJB de 4 variables requiere N=2."))
    steps>=1 || throw(ArgumentError("steps debe ser positivo."))
    a=Float64.(vaxis);b=Float64.(waxis);qs=Float64.(actions)
    length(a)>=2 && length(b)>=2 && length(qs)>=2 ||
        throw(ArgumentError("Grillas insuficientes."))
    all(diff(a).>0) && all(diff(b).>0) ||
        throw(ArgumentError("Los ejes deben crecer estrictamente."))
    all(abs(q)<=s.Umax+1e-12 for q in qs) ||
        throw(ArgumentError("Corrientes fuera de las cotas."))
    axes=(a,b,a,b)
    dims=(length(a),length(b),length(a),length(b))
    V=[fill(Inf,dims) for _ in 0:steps]
    term=V[end]
    for I in CartesianIndices(term)
        x=[axes[j][I[j]] for j in 1:4]
        r2=phase_order_squared(x,s.center)
        if isfinite(r2)
            term[I]=s.kappa*r2
        end
    end
    h=s.T/steps
    for n in steps:-1:1
        W=V[n];next=V[n+1]
        for I in CartesianIndices(W)
            x=[axes[j][I[j]] for j in 1:4]
            if isfinite(phase_order_squared(x,s.center))
                W[I]=minimum(_small_bellman_action(x,q,s,h,axes,next) for q in qs)
            end
        end
    end
    return SmallDesyncHJB(s,axes,qs,h,V)
end

"""Control discreto q,-q a partir de Bellman en el instante de malla `k`."""
function small_desync_feedback(sol::SmallDesyncHJB,x,k::Integer)
    1<=k<length(sol.values) || throw(BoundsError(sol.values,k+1))
    vals=[_small_bellman_action(x,q,sol.spec,sol.h,sol.axes,sol.values[k+1])
          for q in sol.actions]
    j=argmin(vals)
    if !isfinite(vals[j])
        # Ninguna acción tiene valor de continuación válido en la malla.
        # Elegir actions[1] aquí NO sería una decisión de Bellman, sino
        # un artefacto del orden de la lista (habitualmente q=-Umax).
        # El control neutro se utiliza SOLO como fallback de diagnóstico.
        return (u=[0.0,0.0],covered=false)
    end
    return (u=[sol.actions[j],-sol.actions[j]],covered=true)
end

"""Rollout discreto en el mismo modelo y malla temporal HJB.
Si no se cubre un estado, no se presenta como solución HJB válida.
"""
function rollout_small_desync(sol::SmallDesyncHJB,x0=sol.spec.x0;
        mode::Symbol=:hjb)
    mode in (:hjb,:zero) || throw(ArgumentError("mode=:hjb o :zero."))
    length(x0)==4 || throw(DimensionMismatch("Estado de 4 dimensiones."))
    x=Float64.(x0);states=[copy(x)];U=zeros(2,length(sol.values)-1)
    J=0.0;covered=true;traces=[0.0];ts=collect(range(0.0,sol.spec.T,length=length(sol.values)))
    for k in 1:size(U,2)
        control=mode==:hjb ? small_desync_feedback(sol,x,k) : (u=[0.0,0.0],covered=true)
        covered &= control.covered
        u=control.u;q=u[1]
        U[:,k] .= u
        r2=phase_order_squared(x,sol.spec.center)
        isfinite(r2) || throw(DomainError(x,"La fase no está definida."))
        J+=sol.h*(sol.spec.alpha_d*r2+sol.spec.gamma_d*sum(abs2,u))
        x=_small_step(x,q,sol.spec,sol.h)
        push!(states,copy(x));push!(traces,J)
    end
    final_r2=phase_order_squared(x,sol.spec.center)
    isfinite(final_r2) || throw(DomainError(x,"Fase terminal indefinida."))
    J+=sol.spec.kappa*final_r2
    return (t=ts,states=reduce(hcat,states),currents=U,
            r2=[phase_order_squared(y,sol.spec.center) for y in states],
            E=[synchronization_error(y) for y in states],
            running_J=traces,final_J=J,covered=covered,terminal_r2=final_r2)
end
