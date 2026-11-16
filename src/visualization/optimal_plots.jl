# Visualizaciones Etapa E: consumen resultados ya calculados, sin reintegrar EDO.

"""Gráfico comparativo de sincronización: E, r, lambda y costo acumulado."""
function plot_optimal_sync(result,s::SyncOCPSpec,sim,cost;
        comparison=nothing,savepath=nothing)
    fig=Figure(size=(1250,940))
    ax1=Axis(fig[1,1],title="Error de sincronización",xlabel="t",ylabel="E(t)")
    ax2=Axis(fig[1,2],title="Orden geométrico",xlabel="t",ylabel="r(t)")
    ax3=Axis(fig[2,1],title="Acoplamiento por tramos",xlabel="t",ylabel="λ(t)")
    ax4=Axis(fig[2,2],title="Costo acumulado J",xlabel="t",ylabel="∫ L_s dt")
    m=sim.metrics
    lines!(ax1,m.t,m.E;linewidth=2.6,color=:royalblue,label="candidato")
    hlines!(ax1,[s.delta];linestyle=:dash,color=:gray50,label="δ")
    lines!(ax2,m.t,m.r;linewidth=2.6,color=:darkgreen)
    stairs!(ax3,collect(range(0,result.T,length=result.M+1)),
            vcat(result.controls,last(result.controls));step=:post,
            linewidth=2.6,color=:darkorange)
    hlines!(ax3,[s.lambda_ref];color=:gray50,linestyle=:dash)
    lines!(ax4,cost.t,cost.accumulated_cost;linewidth=2.6,color=:purple)
    if comparison !== nothing
        # NamedTuple (metrics=...,label=...), desde simulación independiente.
        lines!(ax1,comparison.metrics.t,comparison.metrics.E;
               label=comparison.label,linestyle=:dash,linewidth=2,color=:tomato)
        axislegend(ax1;position=:rt)
    end
    ylims!(ax2,-0.05,1.05)
    return _save_figure(fig,savepath)
end

"""Mostrar corrientes óptimas candidatas y métricas asociadas."""
function plot_optimal_desync(result,s::DesyncOCPSpec,sim,cost;
        savepath=nothing)
    fig=Figure(size=(1300,940))
    ax1=Axis(fig[1,1],title="Corrientes diferenciales",xlabel="t",ylabel="u_i(t)")
    ax2=Axis(fig[1,2],title="Fases geométricas",xlabel="t",ylabel="r²(t)")
    ax3=Axis(fig[2,1],title="Error de estados",xlabel="t",ylabel="E(t)")
    ax4=Axis(fig[2,2],title="Costo acumulado",xlabel="t",ylabel="∫ L_d dt")
    ts=collect(range(0,s.T,length=result.M+1))
    for i in 1:s.N
        stairs!(ax1,ts,vcat(vec(result.controls[i,:]),result.controls[i,end]);
                step=:post,linewidth=1.6)
    end
    hlines!(ax1,[-s.Umax,s.Umax];linestyle=:dash,color=(:gray50,0.6))
    m=sim.metrics
    lines!(ax2,m.t,m.r2;linewidth=2.6,color=:darkorange)
    lines!(ax3,m.t,m.E;linewidth=2.6,color=:royalblue)
    lines!(ax4,cost.t,cost.accumulated_cost;linewidth=2.6,color=:purple)
    ylims!(ax2,-0.05,1.05)
    return _save_figure(fig,savepath)
end

"""Comparación control aplicado vs ley PMP reconstruida desde los adjuntos."""
function plot_pmp_sync(diag;savepath=nothing)
    fig=Figure(size=(1100,760))
    ax1=Axis(fig[1,1],title="Acoplamiento utilizado y minimizador PMP",
             xlabel="t",ylabel="λ")
    ax2=Axis(fig[2,1],title="Hamiltoniano de sincronización",xlabel="t",ylabel="H_s")
    lines!(ax1,diag.t,diag.lambda;label="λ candidato",linewidth=2.5,color=:royalblue)
    lines!(ax1,diag.t,diag.lambda_pmp;label="λ proyectado PMP",linewidth=2,
           linestyle=:dash,color=:tomato)
    lines!(ax2,diag.t,diag.H;linewidth=2,color=:purple)
    hlines!(ax2,[0.0];linestyle=:dash,color=:gray50)
    axislegend(ax1;position=:rt)
    return _save_figure(fig,savepath)
end

function plot_pmp_desync(diag;savepath=nothing,max_neurons::Integer=6)
    fig=Figure(size=(1140,790))
    ax1=Axis(fig[1,1],title="Corrientes aplicadas (continuo) vs PMP (discontinuo)",
             xlabel="t",ylabel="u_i")
    ax2=Axis(fig[2,1],title="Hamiltoniano de desincronización",xlabel="t",ylabel="H_d")
    for i in 1:min(max_neurons,size(diag.U,1))
        lines!(ax1,diag.t,vec(diag.U[i,:]);linewidth=1.7)
        lines!(ax1,diag.t,vec(diag.U_pmp[i,:]);linewidth=1.2,linestyle=:dash)
    end
    lines!(ax2,diag.t,diag.H;linewidth=2.3,color=:purple)
    return _save_figure(fig,savepath)
end
