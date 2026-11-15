# Ejecutar en la raíz del repositorio:
#   julia --project=. scripts/run_stage_B.jl
# Requiere JLD2: julia --project=. -e 'using Pkg; Pkg.add("JLD2")'

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl
using JLD2
using Dates
using LinearAlgebra

root = normpath(joinpath(@__DIR__, ".."))
figdir = joinpath(root, "docs", "figures")
datadir = joinpath(root, "data", "processed")
mkpath(figdir)
mkpath(datadir)

p = FHNParameters()  # Mismo ejemplo configurable que Etapa A.
h = hopf_analysis(p)
println("Parámetros: ", p)
println("Hopf: I- = ", h.I_minus, ", I+ = ", h.I_plus)
println("Coeficiente radial: κ = ", h.kappa, " (", h.kind, ")")
println("Frecuencia angular crítica: ", h.omega,
        "; período lineal: ", h.period_linear)

# 1) El equilibrio y su estabilidad se pueden recalcular sin integrar EDO.
eq_sweep = equilibrium_sweep(p, range(-0.3, 3.1, length=220))
plot_hopf_diagram(p, eq_sweep;
    savepath=joinpath(figdir, "B_01_hopf_equilibrios.png"))

# 2) Detección desde trayectoria, descartando transitorio.
settings = (transient=400.0, window=350.0, dt=0.25, reltol=1e-9,
            abstol=1e-9, u0=(-1.4,-0.2))
det = detect_limit_cycle(p; settings...)
println("Detección de órbita: ", det.status)
println("Cruces: ", det.crossings, "; T estimado: ", det.period)

det.status == :periodic_candidate ||
    error("No se detectó un candidato periódico. Aumenta transient/window y reintenta.")

orb = refine_periodic_orbit(p, det)
orb.status == :converged ||
    error("Shooting no convergió: $(orb.status), residuo $(orb.residual). " *
          "Prueba más tolerancia temporal o ajusta la conjetura inicial.")
println("Shooting: T = ", orb.period, ", amplitud = ", orb.amplitude_v)
println("Residuo: ", orb.residual, "; Floquet = ", orb.floquet)
println("Estabilidad orbital: ", orb.stability)
plot_periodic_orbit(p, orb;
    savepath=joinpath(figdir, "B_02_ciclo_refinado.png"))

# 3) Estudio de corrientes: NO convierte ausencias de detección en amplitud cero.
scan_currents = collect(range(-0.25, 3.05, length=17))
scan = [detect_limit_cycle(with_current(p, I);
             transient=400.0, window=350.0, dt=0.35)
        for I in scan_currents]
plot_cycle_scan(scan;
    savepath=joinpath(figdir, "B_03_amplitud_periodo.png"))
println("Candidatos de ciclo en barrido: ",
        count(d->d.status == :periodic_candidate, scan), "/", length(scan))

# 4) Continuación por shooting (rama local). Puede fallar: se reporta cada intento.
branch = continue_periodic_branch(p, orb, [1.35, 1.30, 1.25, 1.20]; maxiter=18)
branch_with_seed = [(I=p.I, orbit=orb); branch]
plot_continuation(branch_with_seed;
    savepath=joinpath(figdir, "B_04_continuacion_ciclos.png"))
println("Rama: ", [(r.I, r.orbit.status) for r in branch_with_seed])

# Datos planos, portables entre versiones de este código, sin serializar ODESolution.
outfile = joinpath(datadir, "etapa_B_resultados.jld2")
JLD2.jldsave(outfile;
    generated_at=string(now()),
    parameters=(a=p.a, b=p.b, tau=p.tau, I=p.I),
    integrator="OrdinaryDiffEq.Tsit5",
    initial_condition=collect(settings.u0),
    integration_settings=(transient=settings.transient, window=settings.window,
                          dt=settings.dt, reltol=settings.reltol, abstol=settings.abstol),
    hopf=(I_minus=h.I_minus, I_plus=h.I_plus, kappa=h.kappa,
          omega=h.omega, period_linear=h.period_linear, kind=string(h.kind)),
    orbit=(period=orb.period, amplitude_v=orb.amplitude_v,
           amplitude_w=orb.amplitude_w, residual=orb.residual,
           stability=string(orb.stability),
           floquet=orb.floquet, state0=orb.initial_state),
    orbit_times=Float64.(orb.solution.t),
    orbit_states=hcat(orb.solution.u...),
    scan_currents=scan_currents,
    scan_statuses=[string(d.status) for d in scan],
    scan_periods=[d.period for d in scan],
    scan_amplitudes=[d.amplitude_v for d in scan],
    branch_currents=[x.I for x in branch_with_seed],
    branch_statuses=[string(x.orbit.status) for x in branch_with_seed],
    branch_periods=[x.orbit.period for x in branch_with_seed],
    branch_amplitudes=[x.orbit.amplitude_v for x in branch_with_seed])

println("Figuras: ", figdir)
println("Resultados guardados: ", outfile)
