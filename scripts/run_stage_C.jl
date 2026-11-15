# Uso desde la RAIZ del repositorio:
#   julia --project=. scripts/run_stage_C.jl
#   julia --project=. scripts/run_stage_C.jl --video
# La opción --video produce docs/videos/C_01_red_sincronizacion.mp4.
# Reutiliza datos de Etapa B para inicializar osciladores DESFASADOS si existen.

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl
using JLD2

root = normpath(joinpath(@__DIR__, ".."))
figdir = joinpath(root, "docs", "figures")
viddir = joinpath(root, "docs", "videos")
datadir = joinpath(root, "data", "processed")
foreach(mkpath, (figdir, viddir, datadir))

p = FHNParameters(a=0.7, b=0.5, tau=12.5, I=1.4)
N = 12
A = complete_adjacency(N)
tspan = (0.0, 240.0)
saveat = 0.25
reltol, abstol = 1e-9, 1e-9
center = equilibrium(p)  # Centro geométrico fijo; debe evitarse z_i = center.

# Una órbita de la Etapa B proporciona fases iniciales físicamente coherentes.
# Si no existe la caché, elegimos datos exploratorios y lo indicamos.
bfile = joinpath(datadir, "etapa_B_resultados.jld2")
if isfile(bfile)
    data = JLD2.load(bfile)
    orbit_states = data["orbit_states"]
    oldpars = data["parameters"]
    if isapprox(oldpars.I, p.I; atol=1e-12) &&
       isapprox(oldpars.a, p.a; atol=1e-12) &&
       isapprox(oldpars.b, p.b; atol=1e-12) &&
       isapprox(oldpars.tau, p.tau; atol=1e-12)
        x0 = zeros(2N)
        for i in 1:N
            k = 1 + round(Int, (i-1)*(size(orbit_states,2)-1)/N)
            x0[2i-1:2i] = orbit_states[:,k]
        end
        initial_method = "fases distribuidas de órbita refinada (caché Etapa B)"
    else
        x0 = distributed_initial_state(N,p)
        initial_method = "distribuidos deterministas (caché B incompatible)"
    end
else
    x0 = distributed_initial_state(N,p)
    initial_method = "distribuidos deterministas (sin caché B)"
end
println("Estado inicial: ", initial_method)
println("Parámetros: ", p, "; N=", N)

# El acoplamiento Aij=1 conserva el signo del enunciado:
# lambda < 0 es atractivo; positivo puede desestabilizar diferencias.
lambdas = [-1.5, -0.3, 0.45]
experiments = NamedTuple[]
for lambda in lambdas
    net = FHNNetwork(N,p; adjacency=A, lambda=lambda)
    sol = simulate_network(net; u0=x0, tspan=tspan, saveat=saveat,
                           reltol=reltol, abstol=abstol)
    metrics = network_metrics(sol, net; center=center)
    label = "λ=$(lambda)"
    push!(experiments, (label=label, lambda=lambda,
                        net=net, sol=sol, metrics=metrics))
    filename = joinpath(datadir, "etapa_C_lambda_$(replace(string(lambda), "."=>"p", "-"=>"m")).jld2")
    save_network_result(filename, sol, net, metrics;
        solver="Tsit5", reltol=reltol, abstol=abstol,
        description="Etapa C / $label / $initial_method")
    println(label, ": E(0)=", round(first(metrics.E), digits=5),
            ", E(T)=", round(last(metrics.E), digits=8),
            ", r(T)=", round(last(metrics.r), digits=5),
            ", fases definidas=", all(metrics.phase_defined))
end

plot_coupling_comparison(experiments;
    savepath=joinpath(figdir, "C_01_comparacion_acoplamientos.png"))
main_case = experiments[1]  # lambda=-1.5, sincronización atractiva
plot_network_timeseries(main_case.sol, main_case.net;
    savepath=joinpath(figdir, "C_02_series_temporales.png"))
plot_network_phase(main_case.sol, main_case.net;
    savepath=joinpath(figdir, "C_03_retratos_fase.png"))
plot_network_metrics(main_case.metrics;
    savepath=joinpath(figdir, "C_04_metricas.png"))
plot_network_snapshot(main_case.sol.u[1], main_case.net;
    title="Instante inicial / neuronas desfasadas",
    savepath=joinpath(figdir, "C_05_grafo_inicial.png"))
plot_network_snapshot(main_case.sol.u[end], main_case.net;
    title="Instante final / acoplamiento atractivo",
    savepath=joinpath(figdir, "C_06_grafo_final.png"))

println("Figuras guardadas en: ", figdir)
println("Datos guardados en: ", datadir)

if "--video" in ARGS
    file = joinpath(viddir, "C_01_red_sincronizacion.mp4")
    animate_network_dashboard(main_case.sol, main_case.net, main_case.metrics,
                              file; nframes=150, framerate=24)
    println("Video generado: ", file)
else
    println("Para crear MP4 ejecuta: julia --project=. scripts/run_stage_C.jl --video")
end
