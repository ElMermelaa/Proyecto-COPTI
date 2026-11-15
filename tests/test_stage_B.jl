using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa B — Hopf, ciclos y shooting" begin
    p = FHNParameters()

    @testset "Coeficiente de Hopf, frecuencia y transversalidad" begin
        h = hopf_analysis(p)
        @test h !== nothing
        @test h.kind == :supercritical
        @test h.kappa ≈ -0.25/(8*(12.5-0.25)) atol=1e-13
        @test h.omega ≈ 0.28 atol=1e-12
        @test h.period_linear ≈ 2pi/0.28 atol=1e-10
        @test h.slope_minus > 0
        @test h.slope_plus < 0
        @test abs(stability_info(with_current(p,h.I_minus)).trace) < 1e-9
        @test abs(stability_info(with_current(p,h.I_plus)).trace) < 1e-9
        @test hopf_analysis(FHNParameters(b=0.0)) === nothing
        @test hopf_analysis(FHNParameters(b=0.5,tau=0.25)) === nothing

        sub = hopf_analysis(FHNParameters(b=0.8,tau=4.0))
        @test sub.kind == :subcritical
        @test sub.kappa > 0
        # Degeneración: (2b-1)tau = b^2 => tau=0.64/0.6 para b=0.8
        deg = hopf_analysis(FHNParameters(b=0.8,tau=0.8^2/(2*0.8-1)))
        @test abs(deg.kappa) < 1e-11
        @test deg.kind == :degenerate
    end

    @testset "Ciclo atractivo: no confundir con un transitorio" begin
        det = detect_limit_cycle(p; transient=300, window=300, dt=0.35)
        @test det.status == :periodic_candidate
        @test det.crossings >= 5
        @test 30 < det.period < 40
        @test det.amplitude_v > 2.0
        @test det.relative_period_spread < 0.02
        @test det.relative_amplitude_spread < 0.05

        cycle = refine_periodic_orbit(p, det)
        @test cycle.status == :converged
        @test cycle.residual < 1e-7
        @test cycle.stability == :stable
        @test abs(cycle.phase_multiplier - 1) < 1e-3
        @test abs(cycle.transverse_multiplier) < 1
        @test abs(cycle.period-det.period) < 0.15
        @test cycle.amplitude_v > 2
        @test norm(cycle.solution.u[end] - cycle.solution.u[1]) < 1e-6

        # Prueba de una etapa de continuación por parámetro.
        records = continue_periodic_branch(p, cycle, [1.39]; maxiter=18)
        @test length(records) == 1
        @test records[1].orbit.status == :converged
        @test records[1].orbit.residual < 1e-7
    end

    @testset "Fuera del intervalo oscilatorio" begin
        q = with_current(p,-0.3)
        det = detect_limit_cycle(q; transient=400, window=230, dt=0.4)
        @test det.status != :periodic_candidate
    end

    @testset "Errores explícitos y predicción normal" begin
        @test_throws ArgumentError hopf_seed(p, :desconocido)
        @test_throws ArgumentError refine_periodic_orbit(p, [0.0], 20.0)
        @test_throws ArgumentError detect_limit_cycle(p; dt=-0.1)
        seed = hopf_seed(p, :minus; offset=0.002)
        @test seed.params.I > hopf_analysis(p).I_minus
        @test seed.amplitude_estimate > 0
    end
end

println("Etapa B: tests completados.")
