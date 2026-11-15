# Etapa C. Estado intercalado x=[v1,w1,v2,w2,...,vN,wN].
# Ecuacion oficial: +(lambda/N)*sum_j A[i,j]*(v_i-v_j).

"""    complete_adjacency(N; self_loops=true)

Matriz todos-con-todos. `self_loops=true` reproduce exactamente A_ij=1
para todos i,j, como en el enunciado. La diagonal no cambia el campo
porque v_i-v_i=0.
"""
function complete_adjacency(N::Integer; self_loops::Bool=true)
    N >= 2 || throw(ArgumentError("Se requieren al menos dos neuronas."))
    A = ones(Float64, N, N)
    if !self_loops
        for i in 1:N
            A[i,i] = 0.0
        end
    end
    return A
end

"""    ring_adjacency(N; neighbors=1)

Conectividad de anillo NO dirigida (extensión exploratoria).
No debe confundirse con la red oficial todos-con-todos.
"""
function ring_adjacency(N::Integer; neighbors::Integer=1)
    N >= 3 || throw(ArgumentError("Un anillo requiere N >= 3."))
    1 <= neighbors <= fld(N, 2) ||
        throw(ArgumentError("neighbors debe pertenecer a 1:floor(N/2)."))
    A = zeros(Float64, N, N)
    for i in 1:N, k in 1:neighbors
        A[i, mod1(i+k, N)] = 1.0
        A[i, mod1(i-k, N)] = 1.0
    end
    return A
end

"""    FHNNetwork(N, neuron=FHNParameters(); adjacency=complete_adjacency(N),
                    lambda=-0.5, currents=nothing)

`lambda` es número real o función de tiempo `lambda(t)`.
`currents` define perturbaciones individuales adicionales a `neuron.I`:
vector de longitud N o función `currents(t)` que devuelve dicho vector.
Por defecto, las N neuronas reciben la misma corriente `neuron.I`.
La matriz A es no negativa (posiblemente no simétrica); las garantías de
campo medio y de Lyapunov del cuaderno valen para A todos-con-todos.
"""
struct FHNNetwork{L,C}
    N::Int
    neuron::FHNParameters
    A::Matrix{Float64}
    lambda::L
    currents::C
end

function FHNNetwork(N::Integer, neuron::FHNParameters=FHNParameters();
                    adjacency=complete_adjacency(N), lambda=-0.5, currents=nothing)
    N >= 2 || throw(ArgumentError("Se requieren N >= 2 neuronas."))
    size(adjacency) == (N,N) ||
        throw(ArgumentError("La matriz de conectividad debe ser de tamaño (N,N)."))
    A = Matrix{Float64}(adjacency)
    all(x -> isfinite(x) && x >= 0, A) ||
        throw(ArgumentError("La matriz debe tener entradas no negativas y finitas."))
    (lambda isa Real || applicable(lambda, 0.0)) ||
        throw(ArgumentError("lambda debe ser un número o función lambda(t)."))
    if lambda isa Real
        isfinite(lambda) || throw(ArgumentError("lambda debe ser finito."))
    end
    u = currents === nothing ? zeros(Float64, N) : currents
    if u isa AbstractVector
        length(u) == N || throw(ArgumentError("currents debe tener N entradas."))
        all(x -> x isa Real && isfinite(x), u) ||
            throw(ArgumentError("Corrientes individuales no finitas."))
        u = Float64.(u)
    elseif !applicable(u, 0.0)
        throw(ArgumentError("currents debe ser un vector o función currents(t)."))
    end
    return FHNNetwork{typeof(lambda),typeof(u)}(Int(N), neuron, A, lambda, u)
end

@inline _law_at(c, t) = c isa Real ? c : c(t)
@inline _current_at(c, t) = c isa AbstractVector ? c : c(t)

"""    network_rhs!(du, x, network, t)

Campo acoplado in-place: costo O(N²) para una matriz arbitraria.
"""
function network_rhs!(du, x, net::FHNNetwork, t)
    N = net.N
    length(x) == 2N && length(du) == 2N ||
        throw(DimensionMismatch("Estado y derivada deben tener 2N componentes."))
    lambda = _law_at(net.lambda, t)
    lambda isa Real && isfinite(lambda) ||
        throw(ArgumentError("La intensidad lambda(t) no es real finita."))
    u = _current_at(net.currents, t)
    u isa AbstractVector && length(u) == N ||
        throw(ArgumentError("currents(t) debe devolver un vector de longitud N."))
    all(z -> z isa Real && isfinite(z), u) ||
        throw(ArgumentError("La corriente individual no es real finita."))
    p = net.neuron
    @inbounds for i in 1:N
        vi = x[2i-1]
        wi = x[2i]
        coupling = zero(vi)
        for j in 1:N
            coupling += net.A[i,j] * (vi - x[2j-1])
        end
        du[2i-1] = fhn_nonlinearity(vi) - wi + p.I + u[i] + lambda*coupling/N
        du[2i] = (vi + p.a - p.b*wi)/p.tau
    end
    return nothing
end

"""    network_rhs(x, net, t=0.0)

Versión sin mutación útil para verificaciones y diferenciación automática.
"""
function network_rhs(x::AbstractVector, net::FHNNetwork, t::Real=0.0)
    du = similar(x)
    network_rhs!(du, x, net, t)
    return du
end

"""    synchronized_state(N, z)

Genera un estado exactamente sincronizado `[v,w,v,w,...]`.
"""
function synchronized_state(N::Integer, z)
    N >= 2 || throw(ArgumentError("N debe ser >= 2."))
    length(z) == 2 || throw(ArgumentError("z debe contener (v,w)."))
    all(q -> q isa Real && isfinite(q), z) ||
        throw(ArgumentError("z debe ser real y finito."))
    return repeat(Float64.(collect(z)), N)
end

"""    distributed_initial_state(N, p; radius_v=1.5, radius_w=0.15,
                                    center=equilibrium(p), phase_shift=0.1)

Estados deterministas dispersos alrededor de `center`. Son condiciones
iniciales exploratorias, NO puntos garantizados de una órbita periódica.
"""
function distributed_initial_state(N::Integer, p::FHNParameters;
        radius_v::Real=1.5, radius_w::Real=0.15,
        center=equilibrium(p), phase_shift::Real=0.1)
    N >= 2 || throw(ArgumentError("N debe ser >= 2."))
    length(center) == 2 || throw(ArgumentError("center debe ser (v,w)."))
    all(isfinite, (radius_v, radius_w, phase_shift, center[1], center[2])) ||
        throw(ArgumentError("Se requieren entradas finitas."))
    radius_v >= 0 && radius_w >= 0 ||
        throw(ArgumentError("Los radios deben ser no negativos."))
    x = zeros(Float64, 2N)
    for i in 1:N
        theta = 2pi*(i-1)/N + phase_shift
        x[2i-1] = center[1] + radius_v*cos(theta)
        x[2i]   = center[2] + radius_w*sin(theta)
    end
    return x
end

"""    initial_state_from_orbit(orbit, N; phases=range(0,1,length=N+1)[1:end-1])

Fases = fracciones de período para construir N estados sobre un ciclo
refinado por shooting en la Etapa B. Más apropiado que puntos arbitrarios
para representar osciladores desfasados.
"""
function initial_state_from_orbit(orbit, N::Integer;
        phases=collect(range(0.0, 1.0, length=N+1))[1:end-1])
    N >= 2 || throw(ArgumentError("N debe ser >= 2."))
    orbit.status == :converged || throw(ArgumentError("Se necesita un ciclo refinado."))
    length(phases) == N || throw(ArgumentError("Se requieren N fases."))
    all(q -> q isa Real && isfinite(q), phases) ||
        throw(ArgumentError("Las fases deben ser números finitos."))
    x = zeros(Float64, 2N)
    for i in 1:N
        phase = mod(Float64(phases[i]), 1.0)
        z = orbit.solution(phase*orbit.period)
        x[2i-1] = z[1]
        x[2i] = z[2]
    end
    return x
end

"""    simulate_network(net; u0, tspan=(0.,200.), saveat=0.2, ...)

Integra las 2N ecuaciones con el mismo `Tsit5()` usado en las etapas A/B.
Devuelve `ODESolution`, sin ejecutar gráficos ni almacenar resultados.
"""
function simulate_network(net::FHNNetwork;
        u0=distributed_initial_state(net.N, net.neuron),
        tspan=(0.0, 200.0), saveat=0.2,
        solver=Tsit5(), reltol::Real=1e-9, abstol::Real=1e-9)
    length(u0) == 2 * net.N ||
        throw(ArgumentError("El estado inicial debe tener 2N componentes."))
    all(z -> z isa Real && isfinite(z), u0) ||
        throw(ArgumentError("El estado inicial debe ser real y finito."))
    length(tspan) == 2 || throw(ArgumentError("tspan debe tener dos entradas."))
    t0, tf = Float64.(tspan)
    isfinite(t0) && isfinite(tf) && tf > t0 ||
        throw(ArgumentError("El intervalo temporal debe ser finito y creciente."))
    saveat isa Real && isfinite(saveat) && saveat > 0 ||
        throw(ArgumentError("saveat debe ser un paso real positivo."))
    isfinite(reltol) && reltol > 0 && isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("Tolerancias inválidas."))
    prob = ODEProblem(network_rhs!, Float64.(collect(u0)), (t0,tf), net)
    return solve(prob, solver; saveat=saveat, reltol=reltol, abstol=abstol)
end

"""    network_state_matrix(sol, N)

Matriz `2N × nt` de una solución previamente calculada.
"""
function network_state_matrix(sol, N::Integer)
    length(sol.u) >= 1 || throw(ArgumentError("La solución está vacía."))
    all(z -> length(z) == 2N, sol.u) ||
        throw(DimensionMismatch("La solución no corresponde a N neuronas."))
    return reduce(hcat, sol.u)
end
