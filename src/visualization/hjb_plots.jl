# Gráficos de funciones de valor aproximadas y políticas discretas.

function plot_hjb_sync(sol::ReducedSyncHJB; d0=[0.65,0.12], savepath=nothing)
    fig=Figure(size=(1150,500))
    a1=Axis(fig[1,1],title="HJB de llegada reducida: W(dᵥ,d𝓌)",xlabel="dᵥ",ylabel="d𝓌")
    vals=map(z -> isfinite(z) ? z : NaN,sol.W)
    hm=heatmap!(a1,sol.vgrid,sol.wgrid,vals,colormap=:viridis)
    θ=range(0,2π,length=200)
    lines!(a1,sqrt(sol.spec.delta).*cos.(θ),sqrt(sol.spec.delta).*sin.(θ),color=:white,linewidth=2.5)
    scatter!(a1,[d0[1]],[d0[2]],color=:red,markersize=12)
    Colorbar(fig[1,2],hm,label="Costo aproximado")
    a2=Axis(fig[1,3],title="Realimentación λ(dᵥ,d𝓌)",xlabel="dᵥ",ylabel="d𝓌")
    hm2=heatmap!(a2,sol.vgrid,sol.wgrid,sol.policy,colormap=:balance,
                 colorrange=(sol.spec.lambda_min,sol.spec.lambda_max))
    lines!(a2,sqrt(sol.spec.delta).*cos.(θ),sqrt(sol.spec.delta).*sin.(θ),color=:black,linewidth=2)
    Colorbar(fig[1,4],hm2,label="λ elegido (discreto)")
    if savepath!==nothing
        mkpath(dirname(string(savepath)));save(string(savepath),fig)
    end
    return fig
end

function plot_hjb_desync(sol::SmallDesyncHJB; savepath=nothing)
    # La HJB es de 4 variables. Esta figura muestra un CORTE 2D y
    # diagnostica explicitamente cuando la malla no cubre V(0,x).
    v,w,_,_=sol.axes
    jw=argmin(abs.(w .- sol.spec.x0[2]))
    at0=sol.values[1][:,jw,:,jw]
    covered=isfinite.(at0)
    nvalid=count(covered)
    ntotal=length(at0)
    atT=sol.values[end][:,jw,:,jw]
    terminal_cost=map(z -> isfinite(z) ? z : NaN,atT)
    dh=rollout_small_desync(sol;mode=:hjb)
    dz=rollout_small_desync(sol;mode=:zero)

    fig=Figure(size=(1260,940))
    Label(fig[0,1:2],
          "Etapa F: diagnostico de HJB desincronizacion N=2 (estado 4D)",
          fontsize=22,tellwidth=false)

    title1=nvalid==0 ?
          "V(0,x) sin cobertura: mapa de celdas validas" :
          "V(0,x) donde Bellman tiene valores finitos ($nvalid/$ntotal)"
    a1=Axis(fig[1,1],title=title1,xlabel="v1",ylabel="v2")
    if nvalid==0
        # Mostrar que NO hay datos validos, en vez de un heatmap vacio.
        heatmap!(a1,v,v,Float64.(covered),colormap=:grays,colorrange=(0.0,1.0))
    else
        heatmap!(a1,v,v,map(z -> isfinite(z) ? z : NaN,at0),colormap=:viridis)
    end
    lines!(a1,v,v,color=:white,linestyle=:dash,linewidth=1.5)
    scatter!(a1,[sol.spec.x0[1]],[sol.spec.x0[3]],
             color=:red,markersize=11)

    a2=Axis(fig[1,2],
        title="Costo terminal V(T,x)=kappa r2 (NO es V(0,x))",
        xlabel="v1",ylabel="v2")
    hmT=heatmap!(a2,v,v,terminal_cost,colormap=:viridis)
    lines!(a2,v,v,color=:white,linestyle=:dash,linewidth=1.5)
    scatter!(a2,[sol.spec.x0[1]],[sol.spec.x0[3]],
             color=:red,markersize=11)
    Colorbar(fig[1,3],hmT,label="Costo terminal")

    a3=Axis(fig[2,1],title="Desincronizacion: r2(t)",
            xlabel="t",ylabel="r2")
    lines!(a3,dh.t,dh.r2,color=:royalblue3,linewidth=2.5,
           label="Bellman (con fallback si falta malla)")
    lines!(a3,dz.t,dz.r2,color=:darkorange2,linewidth=2.5,
           linestyle=:dash,label="u=0")
    ylims!(a3,-0.04,1.04)
    axislegend(a3,position=:lb)

    a4=Axis(fig[2,2],
        title="Accion aplicada: q(t) con (u1,u2)=(q,-q)",
        xlabel="t",ylabel="q")
    stairs!(a4,dh.t[1:end-1],vec(dh.currents[1,:]),
            step=:post,color=:royalblue3,linewidth=2.5)
    hlines!(a4,[0.0],color=:gray50,linestyle=:dash)
    ylims!(a4,-sol.spec.Umax-0.1,sol.spec.Umax+0.1)

    state_label=dh.covered ?
        "Trayectoria dentro de cobertura de Bellman (aproximacion discreta)." :
        "ATENCION: fuera de cobertura; se uso u=0 como fallback, NO como control HJB optimo."
    Label(fig[3,1:2],
          "Cobertura V(0,x) en este corte: $nvalid/$ntotal. "*state_label,
          fontsize=15,tellwidth=false)

    if savepath!==nothing
        mkpath(dirname(string(savepath)));save(string(savepath),fig)
    end
    return fig
end

function plot_hjb_comparison(sh,sc,dh,dz;savepath=nothing)
    fig=Figure(size=(1200,850))
    a=Axis(fig[1,1],title="Sincronización (mismo estado inicial, N=2)",xlabel="t",ylabel="E(t)")
    lines!(a,sc.t,sc.E,label="λ=-1.5",linewidth=2.5,color=:darkorange2)
    lines!(a,sh.t,sh.E,label="Política HJB reducida",linewidth=2.5,color=:royalblue3)
    hlines!(a,[0.03],color=:black,linestyle=:dash,label="δ=0.03")
    axislegend(a,position=:rt)
    b=Axis(fig[1,2],title="Control de acoplamiento",xlabel="t",ylabel="λ(t)")
    if length(sh.t)>1
        stairs!(b,sh.t[1:end-1],sh.lambda,step=:post,label="HJB discreta",color=:royalblue3)
    end
    lines!(b,[0.0,maximum((sh.t[end],sc.t[end]))],[-1.5,-1.5],
           label="Constante",color=:darkorange2)
    axislegend(b)
    c=Axis(fig[2,1],title="Desincronización HJB (N=2, horizonte corto)",xlabel="t",ylabel="r²(t)")
    lines!(c,dz.t,dz.r2,label="u=0",linewidth=2.5,color=:darkorange2)
    lines!(c,dh.t,dh.r2,label="HJB discreta",linewidth=2.5,color=:royalblue3)
    axislegend(c)
    d=Axis(fig[2,2],title="Corrientes HJB balanceadas",xlabel="t",ylabel="uᵢ(t)")
    stairs!(d,dh.t[1:end-1],vec(dh.currents[1,:]),step=:post,label="u₁",color=:royalblue3)
    stairs!(d,dh.t[1:end-1],vec(dh.currents[2,:]),step=:post,label="u₂",color=:darkorange2)
    axislegend(d)
    if savepath!==nothing
        mkpath(dirname(string(savepath)));save(string(savepath),fig)
    end
    return fig
end
