# Gráficos de la etapa B. No resuelven EDO y aceptan resultados precalculados.

"""    plot_hopf_diagram(p, sweep; savepath=nothing)

Muestra la rama de equilibrios, su estabilidad y traza. Las líneas discontinuas
marcan los valores críticos de Hopf; NO son una rama de ciclos.
"""
function plot_hopf_diagram(p::FHNParameters, sweep; savepath=nothing)
    h = hopf_analysis(p)
    fig = Figure(size=(1000, 720))
    ax1 = Axis(fig[1,1], title="Rama de equilibrio y cambios de estabilidad",
               xlabel="Corriente I", ylabel="v*(I)")
    ax2 = Axis(fig[2,1], xlabel="Corriente I", ylabel="tr(J*)")
    inds_stable = findall(==(:stable), sweep.classification)
    inds_unstable = findall(==(:unstable), sweep.classification)
    lines!(ax1, sweep.I, sweep.v; color=(:gray, 0.25), linewidth=1)
    if !isempty(inds_stable)
        scatter!(ax1, sweep.I[inds_stable], sweep.v[inds_stable];
                 color=:royalblue3, markersize=6, label="Equilibrio estable")
    end
    if !isempty(inds_unstable)
        scatter!(ax1, sweep.I[inds_unstable], sweep.v[inds_unstable];
                 color=:orangered3, markersize=6, label="Equilibrio inestable")
    end
    lines!(ax2, sweep.I, sweep.trace; color=:darkcyan, linewidth=2)
    hlines!(ax2, [0.0]; color=:black, linestyle=:dash)
    if h !== nothing
        for val in (h.I_minus, h.I_plus)
            vlines!(ax1, [val]; linestyle=:dash, color=(:black,0.55))
            vlines!(ax2, [val]; linestyle=:dash, color=(:black,0.55))
        end
    end
    axislegend(ax1; position=:lt)
    return _save_figure(fig, savepath)
end

"""    plot_periodic_orbit(p, orbit; savepath=nothing)

Muestra una órbita refinada por shooting y su evolución en un período.
"""
function plot_periodic_orbit(p::FHNParameters, orbit; savepath=nothing)
    orbit.status == :converged ||
        throw(ArgumentError("El resultado no es una órbita periódica refinada"))
    s = orbit.solution
    vs = [u[1] for u in s.u]
    ws = [u[2] for u in s.u]
    zstar = equilibrium(p)
    fig = Figure(size=(1080, 650))
    ax1 = Axis(fig[1:2,1], title="Órbita periódica: I=$(round(p.I,digits=4))",
               xlabel="Potencial v", ylabel="Recuperación w")
    lines!(ax1, vs, ws; linewidth=3, color=:royalblue3)
    scatter!(ax1, [vs[1]], [ws[1]]; marker=:circle,
             markersize=14, color=:orangered3, label="Sección de Poincaré")
    scatter!(ax1, [zstar[1]], [zstar[2]]; marker=:star5,
             color=:black, markersize=19, label="Equilibrio")
    axislegend(ax1; position=:lt)
    ax2 = Axis(fig[1,2], title="Potencial durante un período",
               xlabel="t", ylabel="v(t)")
    lines!(ax2, s.t, vs; color=:royalblue3, linewidth=2.2)
    ax3 = Axis(fig[2,2], xlabel="t", ylabel="w(t)")
    lines!(ax3, s.t, ws; color=:darkorange2, linewidth=2.2)
    Label(fig[3,:],
          "T=$(round(orbit.period,digits=5)), Amplitud v=$(round(orbit.amplitude_v,digits=4)), " *
          "|μ⊥|=$(round(abs(orbit.transverse_multiplier),digits=6)), " *
          "cierre=$(round(orbit.residual,sigdigits=3))";
          fontsize=14)
    return _save_figure(fig, savepath)
end

"""    plot_cycle_scan(scan; savepath=nothing)

Para cada corriente, muestra exclusivamente candidatos detectados y no
convierte fallos de detección en amplitud cero.
"""
function plot_cycle_scan(scan; savepath=nothing)
    fig = Figure(size=(970, 630))
    a1 = Axis(fig[1,1], title="Amplitud de potencial (candidatos periódicos)",
              xlabel="Corriente I", ylabel="max(v) - min(v)")
    a2 = Axis(fig[2,1], title="Períodos estimados por cruces de Poincaré",
              xlabel="Corriente I", ylabel="Período T")
    inds = findall(s -> s.status == :periodic_candidate, scan)
    if !isempty(inds)
        currents = [scan[i].current for i in inds]
        amps = [scan[i].amplitude_v for i in inds]
        periods = [scan[i].period for i in inds]
        scatterlines!(a1, currents, amps; color=:royalblue3, markersize=10)
        scatterlines!(a2, currents, periods; color=:darkorange2, markersize=10)
    end
    return _save_figure(fig, savepath)
end

"""    plot_continuation(records; savepath=nothing)

Visualiza una rama de órbitas CORREGIDAS por shooting; interrumpe las curvas
en cualquier registro que no haya convergido.
"""
function plot_continuation(records; savepath=nothing)
    fig = Figure(size=(980, 670))
    a1 = Axis(fig[1,1], title="Continuación por shooting: amplitud del ciclo",
              xlabel="Corriente I", ylabel="Amplitud de v")
    a2 = Axis(fig[2,1], title="Período de la órbita",
              xlabel="Corriente I", ylabel="T")
    good = [x for x in records if x.orbit.status == :converged]
    if !isempty(good)
        ii = [x.I for x in good]
        aa = [x.orbit.amplitude_v for x in good]
        tt = [x.orbit.period for x in good]
        scatterlines!(a1, ii, aa; color=:royalblue3, markersize=9)
        scatterlines!(a2, ii, tt; color=:darkorange2, markersize=9)
    end
    return _save_figure(fig, savepath)
end
