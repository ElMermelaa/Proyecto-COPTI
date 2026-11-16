# Etapa F. HJB ESTACIONARIA de llegada en una subvariedad EXACTAMENTE invariante
# para N=2 y equilibrio con v*=0. No es la HJB completa de redes N>2.
# z1=(d_v,c_w+d_w), z2=(-d_v,c_w-d_w), con c=(0,c_w).
# E = d_v^2+d_w^2, r^2=0 para d != 0; se detiene antes de d=0 (delta>0).

struct ReducedSyncHJB
    spec::SyncOCPSpec
    vgrid::Vector{Float64}
    wgrid::Vector{Float64}
    actions::Vector{Float64}
    h::Float64
    W::Matrix{Float64}
    policy::Matrix{Float64}
    iterations::Int
    max_bellman_change::Float64
end

"""Reconstruye la red N=2 de la subvariedad invariante."""
function reduced_sync_state(d::AbstractVector, p::FHNParameters)
    length(d) == 2 || throw(DimensionMismatch("d necesita dos componentes."))
    c = equilibrium(p)
    abs(c[1]) <= 1e-9 || throw(ArgumentError("Esta reducción exacta exige equilibrio v*=0."))
    return [d[1],c[2]+d[2],-d[1],c[2]-d[2]]
end

"""Dinámica exacta sobre la subvariedad, con +lambda*(v_i-vbar)."""
function reduced_sync_rhs(d::AbstractVector, lambda::Real, p::FHNParameters)
    length(d) == 2 || throw(DimensionMismatch("d necesita dos componentes."))
    abs(equilibrium(p)[1]) <= 1e-9 || throw(ArgumentError("v*=0 es necesario."))
    return _reduced_rhs_fast(d,lambda,p)
end

# Se usa dentro del solver DESPUÉS de verificar una vez la invariancia.
# Evita resolver el equilibrio en los millones de evaluaciones de Bellman.
function _reduced_rhs_fast(d,lambda,p)
    dv,dw=d
    return [(1+lambda)*dv - dv^3/3 - dw,
            (dv - p.b*dw)/p.tau]
end

_reduced_E(d) = d[1]^2 + d[2]^2

function _reduced_L(d,lambda,s::SyncOCPSpec)
    return s.rho + s.alpha*_reduced_E(d) + s.beta +
           s.gamma*(lambda-s.lambda_ref)^2
end

function _reduced_step(d,lambda,p,h)
    k1=_reduced_rhs_fast(d,lambda,p)
    mid=d .+ (h/2).*k1
    k2=_reduced_rhs_fast(mid,lambda,p)
    return d .+ h.*k2
end

function _sl_interp2(xs,ys,W,x,y)
    (xs[1]<=x<=xs[end] && ys[1]<=y<=ys[end]) || return Inf
    i=clamp(searchsortedlast(xs,x),1,length(xs)-1)
    j=clamp(searchsortedlast(ys,y),1,length(ys)-1)
    tx=clamp((x-xs[i])/(xs[i+1]-xs[i]),0.0,1.0)
    ty=clamp((y-ys[j])/(ys[j+1]-ys[j]),0.0,1.0)
    val=0.0
    for (ii,wx) in ((i,1-tx),(i+1,tx)), (jj,wy) in ((j,1-ty),(j+1,ty))
        a=wx*wy
        if a>1e-14
            isfinite(W[ii,jj]) || return Inf
            val += a*W[ii,jj]
        end
    end
    return val
end

function _reduced_action_value(d,lambda,s,h,xs,ys,W)
    dnext=_reduced_step(d,lambda,s.neuron,h)
    c0=_reduced_E(d)-s.delta
    c1=_reduced_E(dnext)-s.delta
    L=_reduced_L(d,lambda,s)
    if c1 <= 0
        # Si el paso atraviesa la superficie terminal, se integra solo hasta
        # la primera llegada aproximada, sin evaluar la fase en d=0.
        θ=clamp(c0/(c0-c1),0.0,1.0)
        return θ*h*L
    end
    vn=_sl_interp2(xs,ys,W,dnext[1],dnext[2])
    return h*L + vn
end

"""    solve_reduced_sync_hjb(spec; vgrid, wgrid, actions, h, maxiter, tol)

Aproxima `0=min_lambda {Ls + grad(W)*F}` para N=2 restringido a
`z1=c+d, z2=c-d`, con v*=0. Usa iteración semi-lagrangiana y frontera
W=0 cuando E<=delta. Solo admite una malla compacta y controles muestreados.
`Inf` significa punto no resuelto dentro del dominio discreto.
"""
function solve_reduced_sync_hjb(s::SyncOCPSpec;
        vgrid=collect(range(-1.5,1.5,length=51)),
        wgrid=collect(range(-0.8,0.8,length=41)),
        actions=collect(range(s.lambda_min,s.lambda_max,length=9)),
        h=0.18,maxiter=250,tol=1e-6)
    s.N==2 || throw(ArgumentError("La HJB reducida se define para N=2."))
    abs(equilibrium(s.neuron)[1])<=1e-9 ||
        throw(ArgumentError("Se requiere v*=0 para reducción invariante exacta."))
    xs=Float64.(vgrid);ys=Float64.(wgrid);aa=Float64.(actions)
    length(xs)>=3 && length(ys)>=3 && length(aa)>=2 && h>0 && maxiter>=1 ||
        throw(ArgumentError("Malla, acciones o paso inválidos."))
    issorted(xs) && all(diff(xs).>0) && issorted(ys) && all(diff(ys).>0) ||
        throw(ArgumentError("Las grillas deben ser estrictamente crecientes."))
    all(a -> s.lambda_min <= a <= s.lambda_max,aa) ||
        throw(ArgumentError("Acciones fuera de Λ."))
    W=fill(Inf,length(xs),length(ys))
    target=BitMatrix(undef,size(W))
    for i in eachindex(xs),j in eachindex(ys)
        target[i,j]=xs[i]^2+ys[j]^2<=s.delta
        target[i,j] && (W[i,j]=0.0)
    end
    count(target)>0 || throw(ArgumentError("La grilla no contiene puntos del objetivo E<=delta."))
    lastchange=Inf;nit=0
    for it in 1:maxiter
        Wnew=copy(W);greatest=0.0;changed=false
        for i in eachindex(xs),j in eachindex(ys)
            target[i,j] && continue
            d=(xs[i],ys[j])
            # Tupla no sirve para broadcasting, usar vector para RK2.
            dv=[d[1],d[2]]
            candidate=minimum(_reduced_action_value(dv,a,s,h,xs,ys,W) for a in aa)
            if candidate < Wnew[i,j]-1e-12
                old=Wnew[i,j]
                Wnew[i,j]=candidate
                changed=true
                greatest=max(greatest,isfinite(old) ? old-candidate : Inf)
            end
        end
        W=Wnew;nit=it;lastchange=greatest
        !changed && break
        isfinite(greatest) && greatest<tol && break
    end
    policy=fill(NaN,size(W))
    for i in eachindex(xs),j in eachindex(ys)
        target[i,j] && continue
        vals=[_reduced_action_value([xs[i],ys[j]],a,s,h,xs,ys,W) for a in aa]
        k=argmin(vals)
        isfinite(vals[k]) && (policy[i,j]=aa[k])
    end
    return ReducedSyncHJB(s,xs,ys,aa,Float64(h),W,policy,nit,lastchange)
end

function reduced_sync_value(sol::ReducedSyncHJB,d)
    _reduced_E(d)<=sol.spec.delta && return 0.0
    return _sl_interp2(sol.vgrid,sol.wgrid,sol.W,d[1],d[2])
end

"""Control minimizador DISCRETO de Bellman y flag de cobertura de la grilla."""
function reduced_sync_feedback(sol::ReducedSyncHJB,d)
    _reduced_E(d)<=sol.spec.delta && return (lambda=sol.spec.lambda_ref,covered=true)
    vals=[_reduced_action_value(d,a,sol.spec,sol.h,sol.vgrid,sol.wgrid,sol.W)
          for a in sol.actions]
    k=argmin(vals)
    return (lambda=isfinite(vals[k]) ? sol.actions[k] : sol.spec.lambda_min,
            covered=isfinite(vals[k]))
end

"""Evalúa feedback discreta o control constante, terminando al alcanzar E<=δ."""
function rollout_reduced_sync(sol::ReducedSyncHJB,d0;
        mode::Symbol=:hjb,lambda_constant=-1.5,dt=0.08,Tmax=30.0)
    length(d0)==2 && all(isfinite,d0) && dt>0 && Tmax>0 ||
        throw(ArgumentError("d0, dt o Tmax incorrectos."))
    mode in (:hjb,:constant) || throw(ArgumentError("mode=:hjb o :constant."))
    d=Float64.(d0);ts=[0.0];Y=[copy(d)];ls=Float64[];Js=[0.0]
    miss=0;status=:not_reached
    while ts[end]<Tmax-1e-12
        if _reduced_E(d)<=sol.spec.delta
            status=:reached;break
        end
        feedback=mode==:hjb ? reduced_sync_feedback(sol,d) : (lambda=lambda_constant,covered=true)
        !feedback.covered && (miss+=1)
        λ=feedback.lambda;h=min(dt,Tmax-ts[end])
        dn=_reduced_step(d,λ,sol.spec.neuron,h)
        θ=1.0
        if _reduced_E(dn)<=sol.spec.delta
            c0=_reduced_E(d)-sol.spec.delta;c1=_reduced_E(dn)-sol.spec.delta
            θ=clamp(c0/(c0-c1),0.0,1.0)
        end
        hact=h*θ
        d= d .+ θ.*(dn.-d)
        push!(ls,λ);push!(Js,Js[end]+hact*_reduced_L(Y[end],λ,sol.spec))
        push!(ts,ts[end]+hact);push!(Y,copy(d))
        if θ<1.0
            status=:reached;break
        end
    end
    return (t=ts,states=reduce(hcat,Y),lambda=ls,J=Js,status=status,
            E=[_reduced_E(y) for y in Y],uncovered_steps=miss,
            full_states=reduce(hcat,[reduced_sync_state(y,sol.spec.neuron) for y in Y]))
end
