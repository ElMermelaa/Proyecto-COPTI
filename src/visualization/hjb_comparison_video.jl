# Video de comparación de controles de Etapa E: MISMO modelo N=4, misma x0,
# misma escala temporal. No es una comparación HJB-vs-PMP (esta se reserva
# para problemas idénticos después de estudiar la resolución HJB).

function _saved_state(times,states,t)
    t<=times[1] && return Float64.(states[:,1])
    t>=times[end] && return Float64.(states[:,end])
    j=clamp(searchsortedlast(times,t),1,length(times)-1)
    θ=(t-times[j])/(times[j+1]-times[j])
    return (1-θ).*states[:,j].+θ.*states[:,j+1]
end

function _sample_phase(state,center)
    r=phase_order(state,center)
    return isfinite(r) ? r : NaN
end

"""    animate_stage_E_comparison(cache,path; delta=0.03, tmax=14, nframes=120)

Crea un video comparativo de los datos verificados en Etapa E contra λ=-1.5
constante. El control calculado TERMINA en su primer horizonte; su estado
queda congelado visualmente después (NO se continúa extrapolando su control).
El control constante continúa hasta tmax. Necesita cache de Etapa E.
"""
function animate_stage_E_comparison(cache,path;
        delta=0.03,tmax=14.0,nframes=120,framerate=24)
    nframes>=2 && tmax>0 && framerate>0 ||
        throw(ArgumentError("nframes>=2 y tiempos positivos."))
    data=load_optimal_result(cache)
    data["kind"] == "synchronization" ||
        throw(ArgumentError("La caché debe ser de sincronización."))
    pars=data["parameters"]
    p=FHNParameters(a=pars.a,b=pars.b,tau=pars.tau,I=pars.I)
    N=Int(pars.N);x0=Float64.(data["initial_state"])
    center=Float64.(data["center"])
    times=Float64.(data["times"])
    states=Float64.(data["states"])
    length(times)==size(states,2) || throw(DimensionMismatch("Caché incompatible."))
    net=FHNNetwork(N,p;lambda=-1.5)
    ref=simulate_network(net;u0=x0,tspan=(0.0,tmax),saveat=0.05)
    frames=collect(range(0.0,tmax,length=nframes))
    frames[1]>=times[1] || throw(ArgumentError("Tiempo inicial incompatible."))
    cand=reduce(hcat,[_saved_state(times,states,t) for t in frames])
    base=reduce(hcat,[Float64.(ref(t)) for t in frames])
    E1=[synchronization_error(cand[:,k]) for k in 1:nframes]
    E0=[synchronization_error(base[:,k]) for k in 1:nframes]
    r1=[_sample_phase(cand[:,k],center) for k in 1:nframes]
    r0=[_sample_phase(base[:,k],center) for k in 1:nframes]
    # Color = potencial, tamaño = potencial normalizado, misma escala en ambos.
    allv=vcat(vec(cand[1:2:end,:]),vec(base[1:2:end,:]))
    low,high=extrema(allv);d=max(high-low,1e-6)
    sizes(vals)=[16+16*clamp((v-low)/d,0.0,1.0) for v in vals]
    pos=network_layout(N)
    edges=_edge_segments(net.A,pos)
    left=Observable(_node_values(cand[:,1],N))
    right=Observable(_node_values(base[:,1],N))
    lsize=Observable(sizes(left[]));rsize=Observable(sizes(right[]))
    cursor=Observable([0.0]);stamp=Observable("t=0 | variable frente a λ=-1.5")
    fig=Figure(size=(1260,1050))
    Label(fig[0,1:2],stamp,fontsize=24,tellwidth=false)
    ax1=Axis(fig[1,1],title="Control variable (Etapa E)",aspect=DataAspect())
    ax2=Axis(fig[1,2],title="Control constante λ=-1.5",aspect=DataAspect())
    for ax in (ax1,ax2)
        !isempty(edges) && linesegments!(ax,edges,color=(:gray40,0.22),linewidth=1.0)
        xlims!(ax,-1.24,1.24);ylims!(ax,-1.24,1.24)
        hidedecorations!(ax);hidespines!(ax)
    end
    scatter!(ax1,pos;color=left,markersize=lsize,colormap=:plasma,
             colorrange=(low,high),strokecolor=:black,strokewidth=1)
    scatter!(ax2,pos;color=right,markersize=rsize,colormap=:plasma,
             colorrange=(low,high),strokecolor=:black,strokewidth=1)
    aE=Axis(fig[2,1:2],title="Mismo estado inicial | E(t)",xlabel="t",ylabel="E")
    lines!(aE,frames,E1,color=:royalblue3,linewidth=2.5,label="Control variable (congelado tras T)")
    lines!(aE,frames,E0,color=:darkorange2,linewidth=2.5,label="λ=-1.5")
    hlines!(aE,[delta],color=:black,linestyle=:dash)
    vlines!(aE,cursor,color=(:black,.5));axislegend(aE,position=:rt)
    aR=Axis(fig[3,1:2],title="Parámetro de Kuramoto r(t)",xlabel="t",ylabel="r")
    lines!(aR,frames,r1,color=:royalblue3,linewidth=2.5)
    lines!(aR,frames,r0,color=:darkorange2,linewidth=2.5)
    vlines!(aR,cursor,color=(:black,.5));ylims!(aR,-.04,1.04)
    mkpath(dirname(string(path)))
    CairoMakie.record(fig,string(path),1:nframes;framerate=framerate) do k
        x=_node_values(cand[:,k],N);y=_node_values(base[:,k],N)
        left[]=x;right[]=y;lsize[]=sizes(x);rsize[]=sizes(y)
        cursor[]=[frames[k]]
        frozen=frames[k]>times[end]
        stamp[]="t=$(round(frames[k],digits=2)) | variable $(frozen ? "detenido tras T" : "en evolución") / constante activa"
    end
    return abspath(string(path))
end

"""Video paralelo de desincronización Etapa E: u_i(t) calculadas frente a u≡0.
Ambos casos usan el MISMO N, corriente basal, lambda0, x0 y horizonte T.
"""
function animate_stage_E_desync_comparison(cache,path;nframes=140,framerate=24)
    nframes>=2 && framerate>0 || throw(ArgumentError("nframes y framerate inválidos."))
    data=load_optimal_result(cache)
    data["kind"]=="desynchronization" ||
        throw(ArgumentError("Se necesita la caché E de desincronización."))
    pars=data["parameters"]
    p=FHNParameters(a=pars.a,b=pars.b,tau=pars.tau,I=pars.I)
    N=Int(pars.N);x0=Float64.(data["initial_state"])
    center=Float64.(data["center"])
    λ0=Float64(first(data["lambda_samples"]))
    times=Float64.(data["times"]);states=Float64.(data["states"])
    size(states,2)==length(times) || throw(DimensionMismatch("Caché de desincronización incompatible."))
    T=Float64(data["T"])
    net=FHNNetwork(N,p;lambda=λ0)
    ref=simulate_network(net;u0=x0,tspan=(0.0,T),saveat=max(T/250,0.05))
    ts=collect(range(0.0,T,length=nframes))
    udata=reduce(hcat,[_saved_state(times,states,t) for t in ts])
    ubase=reduce(hcat,[Float64.(ref(t)) for t in ts])
    E1=[synchronization_error(udata[:,k]) for k in 1:nframes]
    E0=[synchronization_error(ubase[:,k]) for k in 1:nframes]
    r1=[_sample_phase(udata[:,k],center)^2 for k in 1:nframes]
    r0=[_sample_phase(ubase[:,k],center)^2 for k in 1:nframes]
    nodevals=vcat(vec(udata[1:2:end,:]),vec(ubase[1:2:end,:]))
    vmin,vmax=extrema(nodevals);dv=max(vmax-vmin,1e-8)
    siz(vs)=[16+15*clamp((v-vmin)/dv,0.0,1.0) for v in vs]
    positions=network_layout(N);edges=_edge_segments(net.A,positions)
    va=Observable(_node_values(udata[:,1],N));vb=Observable(_node_values(ubase[:,1],N))
    sa=Observable(siz(va[]));sb=Observable(siz(vb[]))
    marker=Observable([0.0]);header=Observable("t=0 | corrientes diferenciadas frente a u=0")
    fig=Figure(size=(1250,1030))
    Label(fig[0,1:2],header,fontsize=23,tellwidth=false)
    g1=Axis(fig[1,1],title="Corrientes controladas uᵢ(t)",aspect=DataAspect())
    g2=Axis(fig[1,2],title="Corrientes nulas uᵢ=0",aspect=DataAspect())
    for ax in (g1,g2)
        !isempty(edges) && linesegments!(ax,edges;color=(:gray40,.2),linewidth=1.0)
        xlims!(ax,-1.23,1.23);ylims!(ax,-1.23,1.23)
        hidedecorations!(ax);hidespines!(ax)
    end
    scatter!(g1,positions;color=va,markersize=sa,colormap=:plasma,
             colorrange=(vmin,vmax),strokecolor=:black,strokewidth=1)
    scatter!(g2,positions;color=vb,markersize=sb,colormap=:plasma,
             colorrange=(vmin,vmax),strokecolor=:black,strokewidth=1)
    ar=Axis(fig[2,1:2],title="Orden de fase r²(t)",xlabel="t",ylabel="r²")
    lines!(ar,ts,r1,color=:royalblue3,label="Corrientes controladas",linewidth=2.5)
    lines!(ar,ts,r0,color=:darkorange2,label="u=0",linewidth=2.5)
    ylims!(ar,-0.04,1.04);vlines!(ar,marker;color=(:black,.5));axislegend(ar)
    ae=Axis(fig[3,1:2],title="Dispersión E(t) de estados",xlabel="t",ylabel="E")
    lines!(ae,ts,E1,color=:royalblue3,linewidth=2.5)
    lines!(ae,ts,E0,color=:darkorange2,linewidth=2.5)
    vlines!(ae,marker,color=(:black,.5))
    mkpath(dirname(string(path)))
    CairoMakie.record(fig,string(path),1:nframes;framerate=framerate) do k
        v1=_node_values(udata[:,k],N);v2=_node_values(ubase[:,k],N)
        va[]=v1;vb[]=v2;sa[]=siz(v1);sb[]=siz(v2)
        marker[]=[ts[k]]
        header[]="t=$(round(ts[k],digits=2)) | λ=$(round(λ0,digits=2)) en ambas redes"
    end
    return abspath(string(path))
end
