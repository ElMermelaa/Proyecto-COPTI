# Videos científicos a partir de los candidatos de la Etapa E.
# CORRECCIÓN: prolonga las EDOs realmente después de T; NO congela estados.
# Dependencias: JLD2 y CairoMakie (ya utilizadas en el proyecto).
# Ejecutar desde la raíz:
# julia --project=. scripts/run_stage_G_presentation_videos.jl

using JLD2
using CairoMakie
using Statistics
using LinearAlgebra
using Printf

const ROOT_G = normpath(joinpath(@__DIR__, ".."))
const DATA_G = joinpath(ROOT_G, "data", "processed")
const VIDEOS_G = joinpath(ROOT_G, "docs", "videos")
# SEGUNDOS REALES DE LA SIMULACIÓN DESPUÉS DEL CONTROL OPTIMIZADO.
# ¡Estos son los parámetros que debes aumentar para ver más dinámica!
const POST_SYNC_G = 80.0
const POST_DESYNC_G = 80.0

# Duración del MP4 y fluidez. NO cambian cuánto tiempo físico se simula.
const FPS_G = 24
const VIDEO_SECONDS_G = 24

# Control posterior a T: :hold mantiene el último valor del candidato;
# :release retira el control (λ=0 para sync; uᵢ=0 para desync).
const POST_CONTROL_G = :hold

include(joinpath(ROOT_G,"src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl: FHNParameters, FHNNetwork, simulate_network,
    synchronization_error, phase_order

function read_stage_e(path::String, kind::String)
    isfile(path) || error("No existe el archivo $path. Ejecuta primero scripts/run_stage_E.jl")
    db = JLD2.load(path)
    for key in ("states", "times", "center", "error_E", "phase_r2", "kind")
        haskey(db, key) || error("Falta '$key' en $path; se esperaba el formato de la Etapa E")
    end
    db["kind"] == kind || error("El archivo $path no es de tipo '$kind'")
    X = Matrix{Float64}(db["states"])
    times = Vector{Float64}(db["times"])
    center = Vector{Float64}(db["center"])
    E = Vector{Float64}(db["error_E"])
    r = sqrt.(clamp.(Vector{Float64}(db["phase_r2"]), 0.0, 1.0))
    N = size(X,1) ÷ 2
    n = size(X,2)
    size(X,1) == 2N || error("La matriz de estados no tiene 2N filas")
    length(times) == n && length(E) == n && length(r) == n ||
        error("Inconsistencia entre tiempos, estados y métricas en $path")
    length(center) == 2 || error("El centro geométrico debe tener dos componentes")
    issorted(times) || error("Los tiempos no están ordenados")
    if kind == "synchronization"
        haskey(db, "lambda_samples") || error("Falta lambda_samples en $path")
        control = Vector{Float64}(db["lambda_samples"])
        length(control) == n || error("lambda_samples no coincide con times")
    else
        haskey(db, "current_samples") || error("Falta current_samples en $path")
        control = Matrix{Float64}(db["current_samples"])
        size(control) == (N,n) || error("current_samples debe tener dimensión N × número de tiempos")
    end
    return (; X, times, center, E, r, control, N, n)
end

# ---------------- Continuación FÍSICA de la EDO después del horizonte T ----------------
function extend_stage_e(path::String, kind::String;
                        post_time::Real, post_control::Symbol=:hold)
    post_time > 0 || error("post_time debe ser positivo")
    post_control in (:hold, :release) || error("post_control debe ser :hold o :release")
    d = read_stage_e(path, kind)
    db = JLD2.load(path)
    pars = db["parameters"]
    p = FHNParameters(a=pars.a, b=pars.b, tau=pars.tau, I=pars.I)
    Tend = last(d.times)
    xT = Float64.(d.X[:,end])

    if kind == "synchronization"
        # El óptimo candidato sólo está definido hasta T.
        # Después: conservar la última λ, o retirar el acoplamiento.
        lambda_post = post_control == :hold ? Float64(last(d.control)) : 0.0
        net = FHNNetwork(d.N, p; lambda=lambda_post)
        post_u = fill(0.0, d.N)
        println("Continuación sync: λ_después_T=",lambda_post,
                " (modo=",post_control,")")
    else
        # Se mantiene el acoplamiento basal (λ₀=-0.15 en Etapa E).
        λ0 = Float64(first(db["lambda_samples"]))
        post_u = post_control == :hold ? Float64.(d.control[:,end]) : zeros(d.N)
        net = FHNNetwork(d.N, p; lambda=λ0, currents=post_u)
        lambda_post = λ0
        println("Continuación desync: λ₀=",λ0,
                ", u_después_T=",round.(post_u,digits=3),
                " (modo=",post_control,")")
    end

    # Integración real desde el estado final guardado: NO se interpola fuera [0,T].
    sol = simulate_network(net; u0=xT,
                           tspan=(Tend,Tend+Float64(post_time)),
                           saveat=0.2, reltol=1e-9, abstol=1e-9)
    # Guardar puntos igualmente espaciados; el último es t=T+post_time.
    npost = max(2, ceil(Int, post_time/0.2))
    tpost = collect(range(Tend,Tend+Float64(post_time),length=npost+1))[2:end]
    Xpost = reduce(hcat, [Float64.(sol(tt)) for tt in tpost])
    times = vcat(d.times,tpost)
    X = hcat(d.X,Xpost)
    E = [synchronization_error(view(X,:,j)) for j in axes(X,2)]
    r = [phase_order(view(X,:,j),d.center) for j in axes(X,2)]
    control = kind == "synchronization" ?
        vcat(d.control,fill(lambda_post,npost)) :
        hcat(d.control,repeat(reshape(post_u,:,1),1,npost))

    println("Tiempo físico: T=",round(Tend,digits=3),
            " → Tfinal=",round(last(times),digits=3),
            "; E(Tfinal)=",round(last(E),digits=5),
            "; r(Tfinal)=",round(last(r),digits=5))
    return (; X, times, center=d.center, E, r, control,
             N=d.N, n=length(times), T=Tend, post_control)
end

function make_video(path::String, outfile::String;
                    kind::String, post_time::Real, post_control::Symbol=:hold)
    d = extend_stage_e(path, kind; post_time=post_time, post_control=post_control)
    X, t, c, E, r, ctrl, N, nt = d.X, d.times, d.center, d.E, d.r, d.control, d.N, d.n
    Tcontrol = d.T
    v, w = X[1:2:end,:], X[2:2:end,:]
    colors_neurons = [:dodgerblue, :tomato, :seagreen, :purple, :goldenrod, :deeppink, :sienna, :teal]
    neuroncolor(i) = colors_neurons[mod1(i, length(colors_neurons))]

    # Dirección geométrica de cada neurona; el centro es el equilibrio individual.
    qx = fill(NaN, N, nt)
    qy = fill(NaN, N, nt)
    distances = zeros(N, nt)
    for k in 1:nt, i in 1:N
        dx, dy = v[i,k] - c[1], w[i,k] - c[2]
        dist = hypot(dx,dy)
        distances[i,k] = dist
        if dist > 1e-10
            qx[i,k], qy[i,k] = dx/dist, dy/dist
        end
    end
    count(x -> x <= 1e-10, distances) > 0 &&
        @warn "Hay neuronas demasiado próximas al centro: algunas fases geométricas son indefinidas"

    # Layout circular fijo del grafo; las aristas representan todos-con-todos.
    angles = collect(range(0, 2π, length=N+1))[1:end-1]
    gx = 0.85 .* cos.(angles)
    gy = 0.85 .* sin.(angles)
    graph_points = Point2f.(gx, gy)
    vmin, vmax = extrema(v)
    if isapprox(vmin,vmax; atol=1e-12)
        vmax = vmin + 1.0
    end

    fig = Figure(size=(1650, 970), fontsize=18)
    title = kind == "synchronization" ?
        "Sincronización: control variable candidato" :
        "Desincronización: corrientes candidatas"
    Label(fig[0,1:3], title, fontsize=27, font=:bold)

    axgraph = Axis(fig[1,1], title="Grafo neuronal: color = vᵢ(t)",
                   aspect=DataAspect(),
                   xticksvisible=false, yticksvisible=false,
                   xticklabelsvisible=false, yticklabelsvisible=false,
                   xgridvisible=false, ygridvisible=false)
    axphase = Axis(fig[1,2], title="Plano de fase (vᵢ,wᵢ)", xlabel="v", ylabel="w")
    axcircle = Axis(fig[1,3], title="Fases respecto del centro c", xlabel="qₓ", ylabel="qᵧ",
                    aspect=DataAspect())
    axE = Axis(fig[2,1], title="Sincronización de estados", xlabel="tiempo", ylabel="E(t)")
    axr = Axis(fig[2,2], title="Alineamiento geométrico", xlabel="tiempo", ylabel="r(t)")
    axctrl = Axis(fig[2,3], title="Control aplicado", xlabel="tiempo",
                  ylabel=(kind == "synchronization" ? "λ(t)" : "uᵢ(t)"))

    for i in 1:N-1, j in i+1:N
        lines!(axgraph, [gx[i],gx[j]], [gy[i],gy[j]], color=(:gray,0.21), linewidth=1)
    end
    limits!(axgraph,-1.12,1.12,-1.12,1.12)

    # Escala fija en color y tamaño: el parpadeo se interpreta como variación de vᵢ.
    current_colors = Observable(collect(v[:,1]))
    current_sizes = Observable(fill(28.0,N))
    scatter!(axgraph, graph_points;
             color=current_colors, colormap=:balance, colorrange=(vmin,vmax),
             markersize=current_sizes, strokecolor=:black, strokewidth=1)
    for i in 1:N
        text!(axgraph,gx[i]*1.07,gy[i]*1.07;
              text=string(i),align=(:center,:center),fontsize=12)
    end

    for i in 1:N
        lines!(axphase,vec(v[i,:]),vec(w[i,:]),color=(neuroncolor(i),0.20),linewidth=1.5)
    end
    scatter!(axphase,[c[1]],[c[2]],color=:black,marker=:cross,markersize=16)
    current_states = Observable([Point2f(v[i,1],w[i,1]) for i in 1:N])
    scatter!(axphase,current_states,color=[neuroncolor(i) for i in 1:N],markersize=16)

    θ = range(0,2π,length=240)
    lines!(axcircle,cos.(θ),sin.(θ),color=:gray,linewidth=2)
    hlines!(axcircle,[0],color=(:gray,0.35))
    vlines!(axcircle,[0],color=(:gray,0.35))
    current_q = Observable([Point2f(qx[i,1],qy[i,1]) for i in 1:N])
    scatter!(axcircle,current_q,color=[neuroncolor(i) for i in 1:N],markersize=19)
    qbar = Observable(Point2f[(0,0),(mean(qx[:,1]),mean(qy[:,1]))])
    lines!(axcircle,qbar,color=:black,linewidth=4)
    limits!(axcircle,-1.13,1.13,-1.13,1.13)

    # Líneas tenues: futuro de la trayectoria. Líneas gruesas: evolución transcurrida.
    lines!(axE,t,E,color=(:steelblue,0.23),linewidth=2)
    lines!(axr,t,r,color=(:darkorange,0.23),linewidth=2)
    now_t = Observable(t[1:1])
    now_E = Observable(E[1:1])
    now_r = Observable(r[1:1])
    lines!(axE,now_t,now_E,color=:steelblue,linewidth=4)
    lines!(axr,now_t,now_r,color=:darkorange,linewidth=4)
    ylims!(axE,0,max(maximum(E)*1.1,0.01))
    ylims!(axr,0,1.05)
    xlims!(axE,first(t),last(t))
    xlims!(axr,first(t),last(t))
    vlines!(axE,[Tcontrol],color=:firebrick,linestyle=:dash,linewidth=2)
    vlines!(axr,[Tcontrol],color=:firebrick,linestyle=:dash,linewidth=2)

    if kind == "synchronization"
        lines!(axctrl,t,ctrl,color=:purple,linewidth=3)
    else
        for i in 1:N
            lines!(axctrl,t,vec(ctrl[i,:]),color=neuroncolor(i),linewidth=2)
        end
    end
    xlims!(axctrl,first(t),last(t))
    vlines!(axctrl,[Tcontrol],color=:firebrick,linestyle=:dash,linewidth=2)

    cursors = [Observable([t[1]]) for _ in 1:3]
    for (ax, cursor) in zip((axE,axr,axctrl),cursors)
        vlines!(ax,cursor,color=(:black,0.65),linewidth=2,linestyle=:dash)
    end

    summary = Observable(@sprintf("t = %.2f | E = %.5f | r = %.5f | T = %.2f",
                                  t[1],E[1],r[1],Tcontrol))
    Label(fig[3,1:3], summary, fontsize=20, tellwidth=false)
    colsize!(fig.layout,1,Relative(1/3))
    colsize!(fig.layout,2,Relative(1/3))
    colsize!(fig.layout,3,Relative(1/3))

    # VIDEO_SECONDS_G segundos de video mostrando TIEMPO FÍSICO real 0→T+post_time.
    # NO se congela el estado ni se reduce el FPS para alargar artificialmente.
    # La interpolación entre muestras sólo suaviza la representación gráfica.
    nframes = FPS_G * VIDEO_SECONDS_G
    nframes >= 2 || error("FPS_G*VIDEO_SECONDS_G debe ser >= 2")
    frame_times = collect(range(first(t),last(t),length=nframes))

    mkpath(dirname(outfile))
    record(fig,outfile,eachindex(frame_times);framerate=FPS_G) do frame
        time_now = frame_times[frame]
        j = clamp(searchsortedlast(t,time_now),1,nt-1)
        h = t[j+1]-t[j]
        theta = h > 0 ? clamp((time_now-t[j])/h,0.0,1.0) : 0.0

        # Interpolación de estados. Se reconstruyen E y r a partir de los
        # estados interpolados para no mostrar métricas incompatibles.
        vt = (1-theta).*v[:,j] .+ theta.*v[:,j+1]
        wt = (1-theta).*w[:,j] .+ theta.*w[:,j+1]
        vbar, wbar = mean(vt), mean(wt)
        Et = sum((vt .- vbar).^2 .+ (wt .- wbar).^2)/N
        qpts = Point2f[]
        qxs = Float64[]
        qys = Float64[]
        for i in 1:N
            dx, dy = vt[i]-c[1], wt[i]-c[2]
            dist = hypot(dx,dy)
            qi = dist > 1e-10 ? Point2f(dx/dist,dy/dist) : Point2f(NaN,NaN)
            push!(qpts,qi)
            push!(qxs,Float64(qi[1]))
            push!(qys,Float64(qi[2]))
        end
        valid = [i for i in 1:N if isfinite(qxs[i]) && isfinite(qys[i])]
        if length(valid) == N
            mx, my = mean(qxs),mean(qys)
            rt = hypot(mx,my)
        else
            mx,my,rt = NaN,NaN,NaN
        end

        current_colors[] = collect(vt)
        current_sizes[] = [20.0+19.0*clamp((vt[i]-vmin)/(vmax-vmin),0,1) for i in 1:N]
        current_states[] = [Point2f(vt[i],wt[i]) for i in 1:N]
        current_q[] = qpts
        qbar[] = Point2f[(0,0),(mx,my)]
        now_t[] = vcat(t[1:j], [time_now])
        now_E[] = vcat(E[1:j], [Et])
        now_r[] = vcat(r[1:j], [rt])
        for cursor in cursors
            cursor[] = [time_now]
        end
        phase_label = time_now <= Tcontrol + 1e-10 ?
            "control candidato" :
            (post_control == :hold ? "continuación: último control fijo" : "continuación: control retirado")
        summary[] = @sprintf("t = %.2f | E = %.5f | r = %.5f | %s | T = %.2f",
                             time_now,Et,rt,phase_label,Tcontrol)
    end
    println("Video generado: ",outfile)
    return nothing
end

function main()
    mkpath(VIDEOS_G)
    sync = joinpath(DATA_G,"etapa_E_sync.jld2")
    desync = joinpath(DATA_G,"etapa_E_desync.jld2")
    post_control = "--release" in ARGS ? :release : POST_CONTROL_G
    mode_tag = post_control == :hold ? "_extendido" : "_libre"
    only_sync = "--sync" in ARGS
    only_desync = "--desync" in ARGS
    if (!only_sync && !only_desync) || only_sync
        make_video(sync,
                   joinpath(VIDEOS_G,"G_01_presentacion_sincronizacion"*mode_tag*".mp4");
                   kind="synchronization",post_time=POST_SYNC_G,post_control=post_control)
    end
    if (!only_sync && !only_desync) || only_desync
        make_video(desync,
                   joinpath(VIDEOS_G,"G_02_presentacion_desincronizacion"*mode_tag*".mp4");
                   kind="desynchronization",post_time=POST_DESYNC_G,post_control=post_control)
    end
    println("Videos con evolución posterior a T en: ",VIDEOS_G)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
