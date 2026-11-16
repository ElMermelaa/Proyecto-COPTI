# Desde la RAÍZ del repositorio:
#   julia --project=. tests/test_stage_D.jl
#   julia --project=. scripts/run_stage_D.jl
#   julia --project=. scripts/run_stage_D.jl --video
# Las figuras y datos se guardan en docs/figures y data/processed.

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl
using JLD2

root = normpath(joinpath(@__DIR__, ".."))
figdir = joinpath(root,"docs","figures")
datadir = joinpath(root,"data","processed")
viddir = joinpath(root,"docs","videos")
foreach(mkpath,(figdir,datadir,viddir))

# Todos estos valores son EXPERIMENTALES, no parámetros fijados por el enunciado.
p = FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
N = 12
A = complete_adjacency(N)
center = equilibrium(p)
tspan = (0.0,240.0)
saveat = 0.25
reltol = 1e-9
abstol = 1e-9
delta = 0.01

# Construimos datos iniciales comunes para hacer comparaciones justas.
bfile = joinpath(datadir,"etapa_B_resultados.jld2")
if isfile(bfile)
    db = JLD2.load(bfile)
    q = db["parameters"]
    fields = (:a,:b,:tau,:I)
    compatible = all(isapprox(getproperty(p,f),getproperty(q,f);atol=1e-12)
                     for f in fields)
    if compatible
        pts = db["orbit_states"]
        x0 = zeros(2N)
        for i in 1:N
            j = 1 + round(Int,(i-1)*(size(pts,2)-1)/N)
            x0[2i-1:2i] .= pts[:,j]
        end
        zsyn = Float64.(pts[:,1])
        println("Condición inicial: fases de ciclo refinado de Etapa B.")
    else
        x0 = distributed_initial_state(N,p)
        zsyn = [-1.4,-0.2]
        println("Caché B incompatible: datos alternativos (no se garantiza ciclo límite inicial).")
    end
else
    x0 = distributed_initial_state(N,p)
    zsyn = [-1.4,-0.2]
    println("Sin caché B: datos exploratorios (no se garantiza ciclo límite inicial).")
end
println("Parámetros: ",p,"; N=",N,"; E(0)=",round(synchronization_error(x0),digits=5))

# ----- Parte 1. Sincronización: acoplamientos constantes y feedback ----
# -1.5 satisface el teorema; -0.15 es ejemplo exploratorio sin garantía global.
# Los pesos J_s son elegidos para estos experimentos y pueden modificarse.
cases = NamedTuple[]
sync_costs = NamedTuple[]
for lambda in (-1.5,-0.3,-0.15,0.45)
    net = FHNNetwork(N,p;adjacency=A,lambda=lambda)
    sol = simulate_network(net;u0=x0,tspan=tspan,saveat=saveat,
                           reltol=reltol,abstol=abstol)
    metrics = network_metrics(sol,net;center=center)
    label = "λ=$(lambda), u=0"
    cost = evaluate_sync_cost(sol,net;delta=delta,rho=1.0,alpha=0.1,
                              beta=0.2,gamma=0.5,lambda_ref=-0.15,center=center)
    push!(cases,(label=label,sol=sol,net=net,metrics=metrics))
    push!(sync_costs,(label=label,cost=cost))
    tag = replace(string(lambda),"-"=>"m","."=>"p")
    save_control_result(joinpath(datadir,"etapa_D_sync_$(tag).jld2"),sol,net,metrics;
         kind="coupling",label=label,sync_cost=cost,reltol=reltol,abstol=abstol)
    println(label,": E(T)=",round(last(metrics.E),digits=7),
            ", r(T)=",round(last(metrics.r),digits=5),
            ", T llegada=",cost.arrival_time,
            ", J_llegada=",cost.arrival_J)
end

# Control λ(t) de una familia cuyo máximo está demostrado (no muestreado):
modulated = AttractiveModulatedCoupling(-1.35,0.25,0.08)
modnet = FHNNetwork(N,p;adjacency=A,lambda=modulated)
modsol = simulate_network(modnet;u0=x0,tspan=tspan,saveat=saveat,
                          reltol=reltol,abstol=abstol)
modmetrics = network_metrics(modsol,modnet;center=center)
modcost = evaluate_sync_cost(modsol,modnet;delta=delta,rho=1.0,alpha=0.1,
                             beta=0.2,gamma=0.5,lambda_ref=-0.15,center=center)
modlabel = "λ(t) ≤ -1.35 (modulado)"
push!(cases,(label=modlabel,sol=modsol,net=modnet,metrics=modmetrics))
push!(sync_costs,(label=modlabel,cost=modcost))
save_control_result(joinpath(datadir,"etapa_D_sync_modulada.jld2"),
    modsol,modnet,modmetrics;kind="bounded_time_coupling",label=modlabel,
    sync_cost=modcost,reltol=reltol,abstol=abstol)
modcert = lyapunov_certificate(modsol,modnet)
println(modlabel,": c=",modcert.c,", E(T)=",last(modmetrics.E),
        ", J_llegada=",modcost.arrival_J)

# Segundo control suficiente: u_i=-k(v_i-v̄) con lambda0=0 y k>1.
base = FHNNetwork(N,p;adjacency=A,lambda=0.0)
feedback = FHNCurrentFeedback(base;k=1.5)
sol_fb = simulate_feedback(feedback;u0=x0,tspan=tspan,saveat=saveat,
                           reltol=reltol,abstol=abstol)
metrics_fb = network_metrics(sol_fb,base;center=center)
append!(cases,[(label="λ₀=0, feedback k=1.5",sol=sol_fb,net=base,
                metrics=metrics_fb)])
energy_fb = total_current_effort(sol_fb,feedback)
save_control_result(joinpath(datadir,"etapa_D_feedback.jld2"),
    sol_fb,feedback,metrics_fb;kind="feedback",label="k=1.5",
    reltol=reltol,abstol=abstol)
println("Feedback: E(T)=",round(last(metrics_fb.E),digits=8),
        ", integral Σuᵢ²=",round(energy_fb,digits=5))

# Certificados teóricos — sólo para casos que satisfacen las HIPÓTESIS.
cert_lambda = lyapunov_certificate(cases[1].sol,cases[1].net)
cert_feedback = lyapunov_certificate(sol_fb,feedback)
println("Certificado λ=-1.5: c=",cert_lambda.c,
        "; exceso máx V-cota=",cert_lambda.max_V_violation)
println("Certificado feedback: c=",cert_feedback.c,
        "; exceso máx V-cota=",cert_feedback.max_V_violation)
println("T suficiente para δ=",delta,": ",
        sufficient_arrival_time(first(cert_lambda.E),delta,p;lambda=-1.5))

plot_lyapunov_certificate(cert_lambda;
    savepath=joinpath(figdir,"D_01_certificado_acoplamiento.png"))
plot_lyapunov_certificate(cert_feedback;
    savepath=joinpath(figdir,"D_02_certificado_feedback.png"))
plot_control_comparison(cases;
    savepath=joinpath(figdir,"D_03_comparacion_sincronizacion.png"))
plot_control_inputs(sol_fb,feedback;
    savepath=joinpath(figdir,"D_04_corrientes_feedback.png"))
plot_sync_cost_traces(sync_costs;
    savepath=joinpath(figdir,"D_05_costos_J_sincronizacion.png"))

# ----- Parte 2. Ruptura de sincronía a partir de un estado idéntico -----
xsync = synchronized_state(N,zsyn)
@assert synchronization_error(xsync) <= 1e-13
# La elección lambda0=-0.15 permite experimentar con corrientes diferenciales
# sin confundir ello con una demostración de desincronización de fase sostenida.
lambda0 = -0.15
profile = balanced_current_profile(N;pattern=:split)
desync_cases = NamedTuple[]
for amp in (0.0,0.25,0.75)
    current = balanced_current_law(profile;amplitude=amp,ramp_time=0.0)
    net = FHNNetwork(N,p;adjacency=A,lambda=lambda0,currents=current)
    sol = simulate_network(net;u0=xsync,tspan=tspan,saveat=saveat,
                           reltol=reltol,abstol=abstol)
    metrics = network_metrics(sol,net;center=center)
    cost = evaluate_desync_cost(sol,net;alpha_d=1.0,gamma_d=0.2,
                                kappa=2.0,Umax=1.0,center=center)
    label = "Amplitud corriente = $amp"
    push!(desync_cases,(label=label,sol=sol,net=net,metrics=metrics,cost=cost))
    tag = replace(string(amp),"."=>"p")
    save_control_result(joinpath(datadir,"etapa_D_desync_$(tag).jld2"),sol,net,metrics;
       kind="desynchronization",label=label,desync_cost=cost,
       reltol=reltol,abstol=abstol)
    println(label,": E(T)=",round(last(metrics.E),digits=6),
            ", r(T)=",round(last(metrics.r),digits=6),
            ", J_desync=",round(cost.total_J,digits=5),
            ", max|uᵢ|=",round(cost.max_current,digits=3))
end

plot_control_comparison(desync_cases;
    savepath=joinpath(figdir,"D_06_ruptura_sincronia.png"))
plot_control_inputs(desync_cases[end].sol,desync_cases[end].net;
    savepath=joinpath(figdir,"D_07_corrientes_desincronizacion.png"))
plot_desync_cost_traces(desync_cases;
    savepath=joinpath(figdir,"D_08_costos_J_desincronizacion.png"))
plot_network_snapshot(desync_cases[end].sol.u[end],desync_cases[end].net;
    title="Diferencias inducidas por corrientes individuales",
    savepath=joinpath(figdir,"D_09_grafo_con_corrientes.png"))

println("Figuras: ",figdir)
println("Datos JLD2: ",datadir)
if "--video" in ARGS
    # La animación recibe la misma solución controlada sin reintegrarla.
    lastcase=desync_cases[end]
    video=joinpath(viddir,"D_01_corrientes_diferenciadas.mp4")
    animate_network_dashboard(lastcase.sol,lastcase.net,lastcase.metrics,video;
                              nframes=160,framerate=24)
    println("Video: ",video)
else
    println("Video opcional: julia --project=. scripts/run_stage_D.jl --video")
end
