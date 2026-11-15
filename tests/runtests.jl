using Test
using LinearAlgebra
using ForwardDiff

# Funciona con el Project.toml ya existente, sin exigir convertirlo en paquete.
include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa A - FitzHugh–Nagumo" begin
    @testset "Validación de parámetros" begin
        p = FHNParameters()
        @test p.a == 0.7
        @test p.b == 0.5
        @test p.tau == 12.5
        @test p.I == 1.4
        @test with_current(p, 2.0).I == 2.0
        @test with_parameter(p, :a, 0.9).a == 0.9
        @test with_parameter(p, :tau, 2.0).tau == 2.0
        @test_throws ArgumentError with_parameter(p, :lambda, 2.0)
        @test_throws ArgumentError FHNParameters(a=-0.1)
        @test_throws ArgumentError FHNParameters(b=1.0)
        @test_throws ArgumentError FHNParameters(tau=0.0)
        @test_throws ArgumentError FHNParameters(I=Inf)
    end

    @testset "Campo vectorial y nulclinas" begin
        p = FHNParameters()
        z = [0.2, 0.4]
        du = zeros(2)
        fhn!(du, z, p, 0.0)
        @test du ≈ fhn_rhs(z, p)
        @test nullcline_v(0.3, p) ≈ 0.3 - 0.3^3/3 + p.I
        @test nullcline_w(0.3, p) ≈ (0.3 + p.a) / p.b
        @test_throws ArgumentError fhn_rhs([0.1], p)
    end

    @testset "Equilibrio único: b>0 y b=0" begin
        p = FHNParameters()  # a=0.7, b=0.5, I=1.4 -> (v*,w*)=(0,1.4)
        z = equilibrium(p)
        @test isapprox(z[1], 0.0; atol=1e-10)
        @test isapprox(z[2], 1.4; atol=1e-10)
        @test equilibrium_residual(p) < 1e-10

        for I in (-2.0, -0.2, 0.4, 1.4, 3.5)
            q = with_current(p, I)
            @test equilibrium_residual(q) < 1e-9
        end

        p0 = FHNParameters(a=0.35, b=0.0, tau=3.0, I=1.1)
        z0 = equilibrium(p0)
        @test z0[1] == -p0.a
        @test z0[2] ≈ fhn_nonlinearity(-p0.a) + p0.I
        @test equilibrium_residual(p0) < 1e-11
        @test_throws DomainError nullcline_w(0.0, p0)
        @test hopf_currents(p0) === nothing
    end

    @testset "Jacobiano y estabilidad" begin
        p = FHNParameters()
        for z in ([0.0, 1.4], [-1.3, 0.1], [0.8, 2.0])
            Jexacto = jacobian_fhn(z[1], p)
            Jad = ForwardDiff.jacobian(u -> fhn_rhs(u, p), z)
            @test Jexacto ≈ Jad atol=1e-11
        end
        info = stability_info(p)
        @test info.classification == :unstable
        @test info.trace > 0
        @test info.determinant > 0
        @test sum(info.eigenvalues) ≈ info.trace
        @test prod(info.eigenvalues) ≈ info.determinant
        @test stability_info(with_current(p, -0.3)).classification == :stable
        @test stability_info(with_current(p, 3.0)).classification == :stable
    end

    @testset "Candidatos de Hopf y barrido" begin
        p = FHNParameters()
        h = hopf_currents(p)
        @test h !== nothing
        @test h.I_minus < h.I_plus
        @test h.omega > 0
        @test abs(stability_info(with_current(p, h.I_minus)).trace) < 1e-9
        @test abs(stability_info(with_current(p, h.I_plus)).trace) < 1e-9
        @test hopf_currents(FHNParameters(tau=0.25)) === nothing

        sw = equilibrium_sweep(p, -0.5:0.5:3.5)
        @test length(sw.I) == length(sw.v) == length(sw.trace)
        @test all(x -> x > 0, sw.determinant)
        for (symbol, vals) in [(:a, [0.3, 0.7, 1.1]),
                               (:b, [0.2, 0.5, 0.8]),
                               (:tau, [5.0, 12.5, 20.0]),
                               (:I, [-0.2, 1.4, 3.0])]
            out = parameter_sweep(p, symbol, vals)
            @test length(out.values) == length(out.v) == length(out.trace)
            @test all(x -> x > 0, out.determinant)
        end
    end

    @testset "Integración y convergencia numérica" begin
        p = FHNParameters()
        zstar = equilibrium(p)
        sol_eq = simulate_individual(p; u0=zstar, tspan=(0.0, 3.0), saveat=0.1)
        @test maximum(norm(u - zstar) for u in sol_eq.u) < 1e-7

        sol1 = simulate_individual(p; u0=(-1.2, 0.2),
                                   tspan=(0.0, 20.0), saveat=0.2,
                                   reltol=1e-7, abstol=1e-7)
        sol2 = simulate_individual(p; u0=(-1.2, 0.2),
                                   tspan=(0.0, 20.0), saveat=0.2,
                                   reltol=1e-10, abstol=1e-10)
        @test norm(sol1.u[end] - sol2.u[end]) < 1e-4
        @test length(sol1.t) > 2
        @test length(simulate_ensemble(p, [[0.0, 1.4], [0.2, 1.2]];
                                       tspan=(0.0, 2.0))) == 2
        @test_throws ArgumentError simulate_individual(p; tspan=(2.0, 1.0))
        @test_throws ArgumentError simulate_individual(p; u0=[1.0])
    end
end

println("Etapa A: pruebas completadas.")
