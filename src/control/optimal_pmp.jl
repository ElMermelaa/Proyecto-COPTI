# MA4703 - Etapa E: Hamiltonianos y condiciones necesarias de Pontryagin.
# No se modifica el signo del acoplamiento del enunciado: +lambda*(v_i - mean(v)).
# Estado y adjunto intercalados: (v1,w1,...,vN,wN).

struct SyncOCPSpec
    neuron::FHNParameters
    N::Int
    x0::Vector{Float64}
    center::Vector{Float64}
    delta::Float64
    lambda_min::Float64
    lambda_max::Float64
    lambda_ref::Float64
    rho::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
end

struct DesyncOCPSpec
    neuron::FHNParameters
    N::Int
    x0::Vector{Float64}
    center::Vector{Float64}
    T::Float64
    lambda0::Float64
    Umax::Float64
    alpha_d::Float64
    gamma_d::Float64
    kappa::Float64
end

function _ocp_data(neuron, N, x0, center)
    N >= 2 || throw(ArgumentError("N debe ser al menos 2."))
    length(x0) == 2N || throw(DimensionMismatch("x0 requiere 2N entradas."))
    length(center) == 2 || throw(DimensionMismatch("center debe tener 2 entradas."))
    all(isfinite,x0) && all(isfinite,center) ||
        throw(ArgumentError("Estado inicial y centro deben ser finitos."))
    return Float64.(collect(x0)), Float64.(collect(center))
end

function SyncOCPSpec(neuron::FHNParameters, N::Integer, x0;
        center=equilibrium(neuron), delta=0.01, lambda_bounds=(-2.5,0.5),
        lambda_ref=-0.15, rho=1.0, alpha=0.1, beta=0.2, gamma=0.5)
    x,c = _ocp_data(neuron,N,x0,center)
    length(lambda_bounds) == 2 || throw(ArgumentError("lambda_bounds=(mín,máx)."))
    lo,hi = Float64.(lambda_bounds)
    all(isfinite,(lo,hi,lambda_ref,rho,alpha,beta,gamma,delta)) ||
        throw(ArgumentError("Todos los parámetros deben ser finitos."))
    lo < hi || throw(ArgumentError("Las cotas de lambda no son válidas."))
    all(q -> q > 0, (rho,alpha,beta,gamma,delta)) ||
        throw(ArgumentError("Los pesos y delta deben ser positivos."))
    E0=synchronization_error(x)
    E0 > delta || throw(ArgumentError("Psin requiere 0<delta<E(x0)."))
    return SyncOCPSpec(neuron,Int(N),x,c,Float64(delta),lo,hi,Float64(lambda_ref),
                       Float64(rho),Float64(alpha),Float64(beta),Float64(gamma))
end

function DesyncOCPSpec(neuron::FHNParameters, N::Integer, x0;
        center=equilibrium(neuron), T=30.0, lambda0=-0.15,
        Umax=1.0, alpha_d=1.0, gamma_d=0.2, kappa=2.0)
    x,c = _ocp_data(neuron,N,x0,center)
    all(isfinite,(T,lambda0,Umax,alpha_d,gamma_d,kappa)) ||
        throw(ArgumentError("Los parámetros deben ser finitos."))
    all(q -> q > 0, (T,Umax,alpha_d,gamma_d,kappa)) ||
        throw(ArgumentError("T, Umax y pesos deben ser positivos."))
    return DesyncOCPSpec(neuron,Int(N),x,c,Float64(T),Float64(lambda0),
                          Float64(Umax),Float64(alpha_d),Float64(gamma_d),Float64(kappa))
end

"""Gradiente exacto del error E(x)=(1/N)sum||z_i-zbar||²."""
function sync_error_gradient(x::AbstractVector)
    N=div(length(x),2)
    N >= 2 && length(x) == 2N || throw(DimensionMismatch("x requiere 2N entradas."))
    vm=sum(x[1:2:end])/N
    wm=sum(x[2:2:end])/N
    g=similar(x)
    for i in 1:N
        g[2i-1]=2*(x[2i-1]-vm)/N
        g[2i]=2*(x[2i]-wm)/N
    end
    return g
end

"""Gradiente exacto de r². Se lanza DomainError si algún z_i coincide con c."""
function phase_order_squared_gradient(x::AbstractVector, center; center_tol=1e-10)
    N=div(length(x),2)
    length(x)==2N && N>=2 || throw(DimensionMismatch("x requiere 2N entradas."))
    length(center)==2 || throw(DimensionMismatch("center debe tener dos componentes."))
    qv=similar(x,N); qw=similar(x,N); d=similar(x,N)
    for i in 1:N
        dv=x[2i-1]-center[1]; dw=x[2i]-center[2]
        di=sqrt(dv*dv+dw*dw)
        di > center_tol || throw(DomainError(di,"La fase geométrica no está definida en el centro."))
        d[i]=di; qv[i]=dv/di; qw[i]=dw/di
    end
    mv=sum(qv)/N; mw=sum(qw)/N
    g=similar(x)
    for i in 1:N
        dotq=qv[i]*mv+qw[i]*mw
        g[2i-1]=2*(mv-qv[i]*dotq)/(N*d[i])
        g[2i]  =2*(mw-qw[i]*dotq)/(N*d[i])
    end
    return g
end

# Mismo E que la etapa C, pero sin validaciones Float64 en cada paso de AD.
function _ocp_E(x)
    N=div(length(x),2)
    vm=sum(x[1:2:end])/N
    wm=sum(x[2:2:end])/N
    return sum((x[2i-1]-vm)^2+(x[2i]-wm)^2 for i in 1:N)/N
end

# Versión genérica, también acepta ForwardDiff.Dual durante la optimización.
function _ocp_r2(x,center; center_tol=1e-10)
    N=div(length(x),2)
    mv=zero(x[1]); mw=zero(x[1])
    for i in 1:N
        dv=x[2i-1]-center[1]; dw=x[2i]-center[2]
        di=sqrt(dv*dv+dw*dw)
        di > center_tol || throw(DomainError(di,"Fase geométrica indefinida (z_i=c)."))
        mv += dv/di
        mw += dw/di
    end
    return (mv*mv+mw*mw)/N^2
end

"""Campo F0(x)+lambda G(x) para la red completa oficial."""
function sync_ocp_field(x, lambda, neuron::FHNParameters)
    N=div(length(x),2)
    length(x)==2N && N>=2 || throw(DimensionMismatch("x requiere 2N entradas."))
    vm=sum(x[1:2:end])/N
    out=similar(x)
    for i in 1:N
        v=x[2i-1]; w=x[2i]
        out[2i-1]=v-v^3/3-w+neuron.I+lambda*(v-vm)
        out[2i]=(v+neuron.a-neuron.b*w)/neuron.tau
    end
    return out
end

function desync_ocp_field(x,u,lambda0,neuron::FHNParameters)
    N=div(length(x),2)
    length(u)==N || throw(DimensionMismatch("Se necesitan N corrientes."))
    out=sync_ocp_field(x,lambda0,neuron)
    for i in 1:N
        out[2i-1] += u[i]
    end
    return out
end

function sync_running_cost(x,lambda,s::SyncOCPSpec)
    r2=_ocp_r2(x,s.center)
    return s.rho+s.alpha*_ocp_E(x)+
           s.beta*(1-r2)^2+s.gamma*(lambda-s.lambda_ref)^2
end

function desync_running_cost(x,u,s::DesyncOCPSpec)
    return s.alpha_d*_ocp_r2(x,s.center)+s.gamma_d*sum(abs2,u)
end

"""Hamiltoniano normal de sincronización: H=L+p'F (convención de mínimo)."""
sync_hamiltonian(x,p,lambda,s::SyncOCPSpec) =
    sync_running_cost(x,lambda,s)+dot(p,sync_ocp_field(x,lambda,s.neuron))

"""Hamiltoniano normal de desincronización: H=L+p'F."""
desync_hamiltonian(x,p,u,s::DesyncOCPSpec) =
    desync_running_cost(x,u,s)+dot(p,desync_ocp_field(x,u,s.lambda0,s.neuron))

"""Ley PMP para acoplamiento, con proyección sobre [lambda_min,lambda_max]."""
function sync_control_pmp(x,p,s::SyncOCPSpec)
    length(x)==2*s.N && length(p)==2*s.N || throw(DimensionMismatch("x y p requieren 2N entradas."))
    vm=sum(x[1:2:end])/s.N
    S=sum(p[2i-1]*(x[2i-1]-vm) for i in 1:s.N)
    return clamp(s.lambda_ref-S/(2*s.gamma),s.lambda_min,s.lambda_max)
end

"""Proyección euclidiana sobre {sum(u)=0, |u_i|<=Umax}."""
function project_balanced_box(y::AbstractVector,Umax::Real)
    N=length(y)
    N>=2 && isfinite(Umax) && Umax>0 || throw(ArgumentError("Se requiere N>=2 y Umax>0."))
    all(isfinite,y) || throw(ArgumentError("El vector a proyectar debe ser finito."))
    # u_i=clamp(y_i-nu,-U,U); la suma es continua y decreciente en nu.
    a=minimum(y)-Umax-1
    b=maximum(y)+Umax+1
    for _ in 1:85
        m=(a+b)/2
        if sum(clamp(z-m,-Umax,Umax) for z in y)>0
            a=m
        else
            b=m
        end
    end
    nu=(a+b)/2
    return [clamp(z-nu,-Umax,Umax) for z in y]
end

"""Ley PMP: resuelve el QP de corrientes con suma cero y saturaciones."""
function desync_control_pmp(p,s::DesyncOCPSpec)
    length(p)==2*s.N || throw(DimensionMismatch("p requiere 2N entradas."))
    return project_balanced_box([-p[2i-1]/(2*s.gamma_d) for i in 1:s.N],s.Umax)
end

"""Sistema adjunto de sincronización: p'=-partial_x H, con lambda fija al derivar."""
function sync_adjoint(x,p,lambda,s::SyncOCPSpec)
    N=s.N
    length(x)==2N && length(p)==2N || throw(DimensionMismatch("x,p requieren 2N entradas."))
    gradE=sync_error_gradient(x)
    gradR=phase_order_squared_gradient(x,s.center)
    r2=_ocp_r2(x,s.center)
    pvbar=sum(p[1:2:end])/N
    dp=similar(x)
    for i in 1:N
        v=x[2i-1]; pv=p[2i-1]; pw=p[2i]
        dp[2i-1]=-s.alpha*gradE[2i-1]+2*s.beta*(1-r2)*gradR[2i-1]-
                  (1-v*v)*pv-lambda*(pv-pvbar)-pw/s.neuron.tau
        dp[2i]  =-s.alpha*gradE[2i]+2*s.beta*(1-r2)*gradR[2i]+pv+
                  s.neuron.b*pw/s.neuron.tau
    end
    return dp
end

"""Sistema adjunto de desincronización con lambda0 constante."""
function desync_adjoint(x,p,u,s::DesyncOCPSpec)
    N=s.N
    length(x)==2N && length(p)==2N && length(u)==N ||
        throw(DimensionMismatch("x,p requieren 2N; u requiere N."))
    gradR=phase_order_squared_gradient(x,s.center)
    pvbar=sum(p[1:2:end])/N
    dp=similar(x)
    for i in 1:N
        v=x[2i-1]; pv=p[2i-1]; pw=p[2i]
        dp[2i-1]=-s.alpha_d*gradR[2i-1]-(1-v*v)*pv-
                  s.lambda0*(pv-pvbar)-pw/s.neuron.tau
        dp[2i]=-s.alpha_d*gradR[2i]+pv+s.neuron.b*pw/s.neuron.tau
    end
    return dp
end
