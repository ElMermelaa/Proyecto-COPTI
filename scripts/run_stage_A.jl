# Desde la raíz del repositorio:
#   julia --project=. scripts/run_stage_A.jl
# Las figuras se guardan en docs/figures/ sin modificar src/ ni notebooks/.

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

root = normpath(joinpath(@__DIR__, ".."))
outdir = joinpath(root, "docs", "figures")
mkpath(outdir)

p = FHNParameters()
println("Parámetros: ", p)
println("Equilibrio: ", equilibrium(p))
println("Estabilidad: ", stability_info(p).classification)
println("Residuo del equilibrio: ", equilibrium_residual(p))
println("Candidatos de Hopf: ", hopf_currents(p))

sol = simulate_individual(p; u0=(-1.4, -0.2), tspan=(0.0, 120.0))
plot_timeseries(sol; savepath=joinpath(outdir, "A_01_series_temporales.png"))

ics = [(-1.4, -0.2), (1.4, 1.9), (0.15, 1.3)]
sols = simulate_ensemble(p, ics; tspan=(0.0, 120.0), saveat=0.12)
plot_phase_portrait(p, sols; savepath=joinpath(outdir, "A_02_retrato_fase.png"))

sweep = equilibrium_sweep(p, range(-0.5, 3.5, length=160))
plot_equilibrium_sweep(p, sweep;
    savepath=joinpath(outdir, "A_03_estabilidad_vs_corriente.png"))

println("Figuras guardadas en: ", outdir)
