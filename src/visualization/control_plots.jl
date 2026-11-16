# ETAPA D — Las figuras consumen resultados YA calculados.

"""    plot_lyapunov_certificate(certificate; savepath=nothing)

V y E observados frente a las cotas exponenciales del cuaderno.
"""
function plot_lyapunov_certificate(cert; savepath=nothing)
    fig = Figure(size=(1050,740))
    ax1 = Axis(fig[1,1],title="Certificado suficiente de Lyapunov (estados)",
               xlabel="Tiempo",ylabel="V(t)")
    ax2 = Axis(fig[2,1],title="Error medio cuadrático y cota suficiente",
               xlabel="Tiempo",ylabel="E(t)")
    lines!(ax1,cert.t,cert.V; color=:royalblue3,linewidth=2.5,label="V calculado")
    lines!(ax1,cert.t,cert.V_upper; color=:black,linewidth=2.0,
           linestyle=:dash,label="V(0) exp(-2ct)")
    lines!(ax2,cert.t,cert.E; color=:darkorange2,linewidth=2.4,label="E calculado")
    lines!(ax2,cert.t,cert.E_upper; color=:black,linewidth=2.0,
           linestyle=:dash,label="Cτ E(0) exp(-2ct)")
    axislegend(ax1; position=:rt)
    axislegend(ax2; position=:rt)
    return _save_figure(fig,savepath)
end

"""    plot_control_comparison(cases; savepath=nothing)

`cases`: colección de `(label, metrics)`; compara E, V y r.
Admite simulaciones con acoplamiento, feedback y corrientes diferenciadas.
"""
function plot_control_comparison(cases; savepath=nothing)
    isempty(cases) && throw(ArgumentError("Se requieren casos."))
    fig = Figure(size=(1120,800))
    axE = Axis(fig[1,1],title="Diferencias entre neuronas",
               xlabel="Tiempo",ylabel="E(t)")
    axV = Axis(fig[2,1],title="Dispersión de Lyapunov",
               xlabel="Tiempo",ylabel="V(t)")
    axR = Axis(fig[3,1],title="Orden geométrico de fase",
               xlabel="Tiempo",ylabel="r(t)")
    for case in cases
        lines!(axE,case.metrics.t,case.metrics.E;linewidth=2.2,label=case.label)
        lines!(axV,case.metrics.t,case.metrics.V;linewidth=2.0,label=case.label)
        lines!(axR,case.metrics.t,case.metrics.r;linewidth=2.2,label=case.label)
    end
    axislegend(axE;position=:rt)
    ylims!(axR,-0.04,1.04)
    return _save_figure(fig,savepath)
end

"""    plot_control_inputs(sol, system; savepath=nothing, max_neurons=12)

Corrientes adicionales u_i(t) y acoplamiento λ(t) del caso indicado.
"""
function plot_control_inputs(sol, system::Union{FHNNetwork,FHNCurrentFeedback};
        savepath=nothing,max_neurons::Integer=12)
    net = system isa FHNNetwork ? system : system.base
    max_neurons >= 1 || throw(ArgumentError("max_neurons debe ser >=1."))
    U = control_current_matrix(sol,system)
    times = Float64.(sol.t)
    lambda = Float64[_law_at(net.lambda,t) for t in times]
    fig = Figure(size=(1050,700))
    ax1 = Axis(fig[1,1],title="Corrientes individuales adicionales",
               xlabel="Tiempo",ylabel="u_i(t)")
    ax2 = Axis(fig[2,1],title="Intensidad de acoplamiento",
               xlabel="Tiempo",ylabel="λ(t)")
    for i in 1:min(size(U,1),max_neurons)
        lines!(ax1,times,vec(U[i,:]);linewidth=1.7)
    end
    hlines!(ax1,[0.0];color=(:black,0.3),linestyle=:dash)
    lines!(ax2,times,lambda;color=:darkcyan,linewidth=2.5)
    return _save_figure(fig,savepath)
end

"""    plot_sync_cost_traces(cases; savepath=nothing)

`cases` contiene `(label, cost)` de `evaluate_sync_cost`.
Muestra J acumulado hasta el horizonte simulado y marca el tiempo de
primera llegada si existe. No interpreta la comparación como optimalidad.
"""
function plot_sync_cost_traces(cases; savepath=nothing)
    isempty(cases) && throw(ArgumentError("Se requieren casos."))
    fig = Figure(size=(1040,570))
    ax = Axis(fig[1,1],title="Costo acumulado para controles prefijados",
              xlabel="Tiempo",ylabel="J acumulado")
    for case in cases
        cost = case.cost
        lines!(ax,cost.t,cost.accumulated_cost;linewidth=2.3,label=case.label)
        if cost.reached
            scatter!(ax,[cost.arrival_time],[cost.arrival_J];markersize=13,
                     marker=:diamond)
        end
    end
    axislegend(ax;position=:lt)
    return _save_figure(fig,savepath)
end

"""    plot_desync_cost_traces(cases; savepath=nothing)

Muestra el costo acumulado integral; el término terminal κr²(T) se informa
separadamente en la tabla/salida y no se suma antes del tiempo final.
"""
function plot_desync_cost_traces(cases; savepath=nothing)
    fig = Figure(size=(1040,570))
    ax = Axis(fig[1,1],title="Desincronización · costo integral acumulado",
              xlabel="Tiempo",ylabel="∫ L_d dt")
    for case in cases
        lines!(ax,case.cost.t,case.cost.accumulated_cost;
               linewidth=2.5,label=case.label)
    end
    axislegend(ax;position=:lt)
    return _save_figure(fig,savepath)
end
