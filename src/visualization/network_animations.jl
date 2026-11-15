# Animación principal: grafo pulsante + plano de fase + E(t) + r(t).
# Usa CairoMakie (instalado en Etapa A) para evitar dependencias adicionales.

"""    animate_network_dashboard(sol, net, metrics, path;
                                  nframes=150, framerate=24, max_edges=350)

Graba un archivo `.mp4` o `.gif` desde resultados YA calculados.
Paneles: (1) grafo de nodos cuya tonalidad y tamaño siguen a vᵢ,
(2) posiciones en el plano (v,w), (3) E(t), (4) r(t).
Los nodos NO tienen tamaños físicos: el pulso es una codificación visual.
`nframes` solo define la densidad del video; no integra de nuevo la EDO.
"""
function animate_network_dashboard(sol, net::FHNNetwork, metrics, path;
        nframes::Integer=150, framerate::Integer=24, max_edges::Integer=350)
    nframes >= 2 || throw(ArgumentError("nframes debe ser >= 2."))
    framerate > 0 || throw(ArgumentError("framerate debe ser > 0."))
    filename = abspath(string(path))
    ext = lowercase(splitext(filename)[2])
    ext in (".mp4", ".gif") ||
        throw(ArgumentError("El archivo de salida debe tener extensión .mp4 o .gif."))
    length(metrics.t) == length(sol.t) ||
        throw(DimensionMismatch("Las métricas deben corresponder a la solución."))
    N = net.N
    states = network_state_matrix(sol, N)
    times = Float64.(sol.t)
    ids = unique(round.(Int, range(1, length(times), length=min(nframes, length(times)))))
    length(ids) >= 2 || throw(ArgumentError("Se necesitan al menos dos tiempos guardados."))

    vs = states[1:2:end, :]
    ws = states[2:2:end, :]
    vmin, vmax = extrema(vs)
    wmin, wmax = extrema(ws)
    dv = max(vmax-vmin, 1e-8)
    dw = max(wmax-wmin, 1e-8)
    vcolor = (vmin-0.02*dv, vmax+0.02*dv)
    sizes_for(values) = [16.0+14.0*clamp((v-vmin)/dv,0.0,1.0) for v in values]

    positions = network_layout(N)
    vals0 = _node_values(states[:, ids[1]], N)
    nodevals = Observable(vals0)
    sizes = Observable(sizes_for(vals0))
    phasepos = Observable(_phase_positions(states[:, ids[1]], N))
    timemarker = Observable([times[ids[1]]])
    timestamp = Observable("t = $(round(times[ids[1]],digits=2))")

    fig = Figure(size=(1250, 900))
    Label(fig[0,1:2], timestamp; fontsize=23, tellwidth=false)
    axG = Axis(fig[1,1], title="Red neuronal | color y pulso = potencial vᵢ",
               aspect=DataAspect())
    axP = Axis(fig[1,2], title="Plano de fase", xlabel="v", ylabel="w")
    axE = Axis(fig[2,1], title="Error de estados E(t)", xlabel="t", ylabel="E")
    axR = Axis(fig[2,2], title="Orden de fase r(t)", xlabel="t", ylabel="r")

    segments = _edge_segments(net.A, positions; max_edges=max_edges)
    if !isempty(segments)
        linesegments!(axG, segments; color=(:gray40,0.18), linewidth=1.0)
    end
    scatter!(axG, positions; color=nodevals, markersize=sizes,
             colormap=:plasma, colorrange=vcolor,
             strokewidth=1.0, strokecolor=:black)
    xlims!(axG,-1.24,1.24)
    ylims!(axG,-1.24,1.24)
    hidedecorations!(axG)
    hidespines!(axG)

    # Las trayectorias completas se dibujan en gris para dar contexto geométrico.
    for i in 1:min(N, 14)
        lines!(axP, vec(vs[i,:]), vec(ws[i,:]);
               linewidth=1.0, color=(:gray50,0.24))
    end
    scatter!(axP, phasepos; color=nodevals, colormap=:plasma,
             colorrange=vcolor, markersize=17,
             strokewidth=0.75, strokecolor=:black)
    xlims!(axP, vmin-0.1*dv, vmax+0.1*dv)
    ylims!(axP, wmin-0.1*dw, wmax+0.1*dw)

    lines!(axE, metrics.t, metrics.E; linewidth=2.3, color=:darkorange2)
    lines!(axR, metrics.t, metrics.r; linewidth=2.3, color=:purple3)
    vlines!(axE, timemarker; color=(:black,0.55), linewidth=2)
    vlines!(axR, timemarker; color=(:black,0.55), linewidth=2)
    ylims!(axR, -0.04, 1.04)

    mkpath(dirname(filename))
    CairoMakie.record(fig, filename, ids; framerate=framerate) do k
        values = _node_values(states[:, k], N)
        nodevals[] = values
        sizes[] = sizes_for(values)
        phasepos[] = _phase_positions(states[:, k], N)
        timemarker[] = [times[k]]
        timestamp[] = "t = $(round(times[k],digits=2)) | λ = $(_law_at(net.lambda, times[k]))"
    end
    return filename
end
