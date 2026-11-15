# Gráficos estáticos. Todas las funciones reciben soluciones ya calculadas.

"""    network_layout(N)

Posiciones deterministas sobre una circunferencia; no dependen del estado.
"""
function network_layout(N::Integer)
    N >= 2 || throw(ArgumentError("N debe ser >= 2."))
    return Point2f[Point2f(cos(2pi*(i-1)/N + pi/2),
                           sin(2pi*(i-1)/N + pi/2)) for i in 1:N]
end

function _edge_segments(A::AbstractMatrix, positions; max_edges::Integer=350)
    N = size(A,1)
    length(positions) == N || throw(DimensionMismatch("Posiciones incorrectas."))
    max_edges >= 0 || throw(ArgumentError("max_edges debe ser no negativo."))
    pairs = Tuple{Int,Int}[]
    for i in 1:N-1, j in i+1:N
        if A[i,j] > 0 || A[j,i] > 0
            push!(pairs, (i,j))
        end
    end
    # Evitar saturación visual para grafos completos grandes.
    if length(pairs) > max_edges
        indices = unique(round.(Int, range(1, length(pairs), length=max_edges)))
        pairs = pairs[indices]
    end
    segments = Point2f[]
    for (i,j) in pairs
        push!(segments, positions[i], positions[j])
    end
    return segments
end

function _node_values(state, N)
    length(state) == 2N || throw(DimensionMismatch("Estado debe tener 2N entradas."))
    return Float64[state[2i-1] for i in 1:N]
end

function _phase_positions(state, N)
    length(state) == 2N || throw(DimensionMismatch("Estado debe tener 2N entradas."))
    return Point2f[Point2f(state[2i-1], state[2i]) for i in 1:N]
end

"""    plot_network_snapshot(state, net; title, savepath)

Vista de la conectividad, con COLOR del nodo igual al potencial vᵢ.
La posición circular es ilustrativa y no representa distancias físicas.
"""
function plot_network_snapshot(state, net::FHNNetwork;
        title="Red neuronal · potencial instantáneo", savepath=nothing,
        max_edges::Integer=350, colorrange=(-2.5, 2.5))
    N = net.N
    vals = _node_values(state, N)
    pos = network_layout(N)
    fig = Figure(size=(750, 630))
    ax = Axis(fig[1,1], title=title, aspect=DataAspect())
    edges = _edge_segments(net.A, pos; max_edges=max_edges)
    if !isempty(edges)
        linesegments!(ax, edges; color=(:gray45, 0.2), linewidth=1.0)
    end
    nodes = scatter!(ax, pos; color=vals, colorrange=colorrange,
                     colormap=:plasma, markersize=24,
                     strokewidth=1.2, strokecolor=:black)
    Colorbar(fig[1,2], nodes; label="Potencial vᵢ")
    xlims!(ax, -1.23, 1.23)
    ylims!(ax, -1.23, 1.23)
    hidedecorations!(ax)
    hidespines!(ax)
    return _save_figure(fig, savepath)
end

"""    plot_network_timeseries(sol, net; max_neurons=14, savepath=nothing)

Muestra vᵢ(t) y wᵢ(t) de un subconjunto representativo de neuronas.
"""
function plot_network_timeseries(sol, net::FHNNetwork;
        max_neurons::Integer=14, savepath=nothing)
    N = net.N
    max_neurons > 0 || throw(ArgumentError("max_neurons debe ser > 0."))
    count = min(N, max_neurons)
    states = network_state_matrix(sol, N)
    fig = Figure(size=(1120, 690))
    ax1 = Axis(fig[1,1], title="Potenciales neuronales", xlabel="t", ylabel="vᵢ(t)")
    ax2 = Axis(fig[2,1], title="Recuperación", xlabel="t", ylabel="wᵢ(t)")
    for i in 1:count
        lines!(ax1, sol.t, vec(states[2i-1,:]); linewidth=1.6)
        lines!(ax2, sol.t, vec(states[2i,:]); linewidth=1.3)
    end
    return _save_figure(fig, savepath)
end

"""    plot_network_phase(sol, net; max_neurons=14, savepath=nothing)

Trayectorias individuales en el plano (v,w) y posiciones finales.
"""
function plot_network_phase(sol, net::FHNNetwork;
        max_neurons::Integer=14, savepath=nothing)
    N = net.N
    count = min(N, max_neurons)
    states = network_state_matrix(sol, N)
    fig = Figure(size=(850, 650))
    ax = Axis(fig[1,1], title="Retrato de fase de la red", xlabel="v", ylabel="w")
    for i in 1:count
        lines!(ax, vec(states[2i-1,:]), vec(states[2i,:]); linewidth=1.4)
    end
    scatter!(ax, _phase_positions(sol.u[end], N);
             color=:black, markersize=8)
    return _save_figure(fig, savepath)
end

"""    plot_network_metrics(metrics; savepath=nothing)

E(t), V(t), r(t). NaN en r no se reemplaza por sincronía perfecta.
"""
function plot_network_metrics(metrics; savepath=nothing)
    fig = Figure(size=(1040, 750))
    axE = Axis(fig[1,1], title="Error de estados E(t)", xlabel="t", ylabel="E")
    axV = Axis(fig[2,1], title="Dispersión de Lyapunov V(t)", xlabel="t", ylabel="V")
    axR = Axis(fig[3,1], title="Orden de fase r(t)", xlabel="t", ylabel="r")
    lines!(axE, metrics.t, metrics.E; linewidth=2.3, color=:darkorange2)
    lines!(axV, metrics.t, metrics.V; linewidth=2.3, color=:royalblue3)
    lines!(axR, metrics.t, metrics.r; linewidth=2.3, color=:purple3)
    ylims!(axR, -0.04, 1.04)
    return _save_figure(fig, savepath)
end

"""    plot_coupling_comparison(experiments; savepath=nothing)

`experiments` es una colección de (label, metrics). Compara E(t) y r(t)
para intensidades de acoplamiento diferentes con el MISMO estado inicial.
"""
function plot_coupling_comparison(experiments; savepath=nothing)
    fig = Figure(size=(1080, 620))
    axE = Axis(fig[1,1], title="Sincronización de estados", xlabel="t", ylabel="E(t)")
    axR = Axis(fig[2,1], title="Alineamiento geométrico", xlabel="t", ylabel="r(t)")
    for case in experiments
        lines!(axE, case.metrics.t, case.metrics.E;
               linewidth=2.3, label=case.label)
        lines!(axR, case.metrics.t, case.metrics.r;
               linewidth=2.3, label=case.label)
    end
    ylims!(axR, -0.04, 1.04)
    axislegend(axE; position=:rt)
    return _save_figure(fig, savepath)
end
