# Las funciones de este archivo devuelven Figure de Makie y no resuelven EDO.
# Esto permite reconstruir una figura sin repetir una simulación costosa.

function _save_figure(fig, savepath)
    if savepath !== nothing
        path = abspath(string(savepath))
        mkpath(dirname(path))
        CairoMakie.save(path, fig)
    end
    return fig
end

"""    plot_timeseries(sol; title="Evolución temporal", savepath=nothing)

Grafica `v(t)` y `w(t)` desde una solución producida por `simulate_individual`.
"""
function plot_timeseries(sol; title="Evolución temporal", savepath=nothing)
    times = sol.t
    vs = [u[1] for u in sol.u]
    ws = [u[2] for u in sol.u]

    fig = Figure(size=(950, 630))
    ax_v = Axis(fig[1, 1], title=title, xlabel="Tiempo t", ylabel="Potencial v(t)")
    ax_w = Axis(fig[2, 1], xlabel="Tiempo t", ylabel="Recuperación w(t)")
    lines!(ax_v, times, vs, color=:royalblue3, linewidth=2.2)
    lines!(ax_w, times, ws, color=:darkorange2, linewidth=2.2)
    return _save_figure(fig, savepath)
end

"""    plot_phase_portrait(p, solutions=[]; vlims=(-2.5,2.5),
                                wlims=(-0.8,2.8), draw_field=true,
                                savepath=nothing)

Retrato de fase con nulclinas, equilibrio, trayectorias y campo de direcciones.
Para `b=0`, muestra correctamente la nulclina vertical `v=-a`.
`solutions` es una colección de soluciones YA calculadas.
"""
function plot_phase_portrait(p::FHNParameters, solutions::AbstractVector=Any[];
        vlims=(-2.5, 2.5),
        wlims=(-0.8, 2.8),
        draw_field::Bool=true,
        savepath=nothing)

    fig = Figure(size=(970, 690))
    ax = Axis(fig[1, 1], title="Retrato de fase FitzHugh–Nagumo",
              xlabel="Potencial v", ylabel="Recuperación w")

    if draw_field
        points = Point2f[]
        directions = Vec2f[]
        for v in range(vlims[1], vlims[2], length=17),
            w in range(wlims[1], wlims[2], length=13)
            dv, dw = fhn_vector_field(v, w, p)
            if hypot(dv, dw) > 1e-10
                push!(points, Point2f(v, w))
                push!(directions, Vec2f(dv, dw))
            end
        end
        arrows2d!(ax, points, directions;
                  normalize=true, lengthscale=0.16, align=:center,
                  color=(:gray, 0.38))
    end

    vgrid = collect(range(vlims[1], vlims[2], length=400))
    lines!(ax, vgrid, nullcline_v.(vgrid, Ref(p));
           linewidth=3, color=:crimson, label="dv/dt = 0")
    if iszero(p.b)
        vlines!(ax, [-p.a]; color=:darkgreen, linewidth=2.8,
                linestyle=:dash, label="dw/dt = 0 (vertical)")
    else
        lines!(ax, vgrid, nullcline_w.(vgrid, Ref(p));
               linewidth=2.8, color=:darkgreen, linestyle=:dash,
               label="dw/dt = 0")
    end

    palette = (:royalblue3, :darkorange2, :purple4, :darkcyan, :magenta3)
    for (i, sol) in enumerate(solutions)
        vs = [u[1] for u in sol.u]
        ws = [u[2] for u in sol.u]
        color = palette[mod1(i, length(palette))]
        lines!(ax, vs, ws; linewidth=2, color=color, label="Trayectoria $i")
        scatter!(ax, [vs[1]], [ws[1]]; color=color, markersize=9)
    end
    z = equilibrium(p)
    scatter!(ax, [z[1]], [z[2]]; color=:black, marker=:star5,
             markersize=20, label="Equilibrio")

    xlims!(ax, vlims...)
    ylims!(ax, wlims...)
    axislegend(ax; position=:rb, framevisible=true)
    return _save_figure(fig, savepath)
end

"""    plot_nullclines(p; kwargs...)

Retrato con campo y nulclinas, sin integrar ninguna trayectoria.
"""
plot_nullclines(p::FHNParameters; kwargs...) = plot_phase_portrait(p; kwargs...)

"""    plot_equilibrium_sweep(p, sweep; savepath=nothing)

Grafica el equilibrio `v*(I)` y la traza del jacobiano al variar la corriente.
Marca las corrientes candidatas de Hopf cuando están definidas.
"""
function plot_equilibrium_sweep(p::FHNParameters, sweep; savepath=nothing)
    fig = Figure(size=(1000, 660))
    ax1 = Axis(fig[1, 1], title="Equilibrio frente a la corriente externa",
               xlabel="Corriente I", ylabel="v*(I)")
    ax2 = Axis(fig[2, 1], xlabel="Corriente I", ylabel="tr J*(I)")

    lines!(ax1, sweep.I, sweep.v; color=:royalblue3, linewidth=2.6)
    hlines!(ax2, [0.0]; color=:gray40, linestyle=:dash)
    lines!(ax2, sweep.I, sweep.trace; color=:darkorange2, linewidth=2.4)

    hopf = hopf_currents(p)
    if hopf !== nothing
        for ih in (hopf.I_minus, hopf.I_plus)
            vlines!(ax1, [ih]; color=(:black, 0.45), linestyle=:dash)
            vlines!(ax2, [ih]; color=(:black, 0.45), linestyle=:dash)
        end
    end
    return _save_figure(fig, savepath)
end
