using Test
using LinearAlgebra
using ForwardDiff
using JLD2

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa C — redes neuronales y métricas" begin
    p = FHNParameters()
    N = 6
    A = complete_adjacency(N)
    @testset "Conectividad y validación" begin
        @test size(A) == (N,N)
        @test all(A .== 1)
        @test all(diag(complete_adjacency(N; self_loops=false)) .== 0)
        @test all(complete_adjacency(N; self_loops=false)[1,2:end] .== 1)
        @test ring_adjacency(6) == transpose(ring_adjacency(6))
        @test all(sum(ring_adjacency(6), dims=2) .== 2)
        @test_throws ArgumentError complete_adjacency(1)
        @test_throws ArgumentError ring_adjacency(2)
        @test_throws ArgumentError ring_adjacency(6; neighbors=5)
        @test_throws ArgumentError FHNNetwork(N,p; adjacency=rand(N-1,N-1))
        @test_throws ArgumentError FHNNetwork(N,p; adjacency=-ones(N,N))
        @test_throws ArgumentError FHNNetwork(N,p; lambda=Inf)
        @test_throws ArgumentError FHNNetwork(N,p; currents=[1.0])
        @test is_mean_field(FHNNetwork(N,p))
        @test !is_mean_field(FHNNetwork(N,p; adjacency=ring_adjacency(N)))
    end

    @testset "Campo medio, signo y corriente externa" begin
        net = FHNNetwork(N,p; lambda=-0.6)
        x = distributed_initial_state(N,p)
        vm = sum(x[1:2:end])/N
        y = network_rhs(x, net)
        for i in 1:N
            vi, wi = x[2i-1:2i]
            @test y[2i-1] ≈ fhn_nonlinearity(vi)-wi+p.I-0.6*(vi-vm) atol=1e-12
            @test y[2i] ≈ (vi+p.a-p.b*wi)/p.tau atol=1e-12
        end
        @test network_rhs(x, FHNNetwork(N,p; lambda=t->-0.6), 11.0) ≈ y
        @test network_rhs(x, FHNNetwork(N,p; lambda=-0.6,
             currents=t->zeros(N)), 2.0) ≈ y
        local_net = FHNNetwork(N,p; lambda=0.0,
                               currents=[0.5; zeros(N-1)])
        local_y = network_rhs(x, local_net)
        baseline = network_rhs(x, FHNNetwork(N,p; lambda=0.0))
        @test local_y[1] - baseline[1] ≈ 0.5 atol=1e-12
        @test local_y[3:end] ≈ baseline[3:end]
        @test norm(network_rhs(x, FHNNetwork(N,p; lambda=0.45))-
                   network_rhs(x, FHNNetwork(N,p; lambda=-0.6))) > 0
        @test_throws DimensionMismatch network_rhs(zeros(3),net)
        # Derivación automática del campo de red para futuras ecuaciones adjuntas.
        J = ForwardDiff.jacobian(z -> network_rhs(z, net), x)
        @test size(J) == (2N, 2N)
        @test J[1,1] ≈ 1-x[1]^2 + (-0.6)*(1-1/N) atol=1e-10
        @test J[1,3] ≈ -(-0.6)/N atol=1e-10
    end

    @testset "Variedad sincronizada y consistencia del modelo individual" begin
        z0 = [-1.4, -0.2]
        x0 = synchronized_state(N,z0)
        for lambda in (-2.5, -0.4, 0.0, 0.8)
            net = FHNNetwork(N,p; lambda=lambda)
            d = network_rhs(x0,net)
            @test d ≈ repeat(fhn_rhs(z0,p),N) atol=1e-12
            sol = simulate_network(net; u0=x0, tspan=(0.,50.), saveat=0.25)
            individual = simulate_individual(p; u0=z0, tspan=(0.,50.), saveat=0.25)
            @test maximum(synchronization_error(z) for z in sol.u) < 1e-16
            @test maximum(norm(sol.u[k][1:2] - individual.u[k]) for k in eachindex(sol.u)) < 2e-5
        end
    end

    @testset "Métricas y cota de Lyapunov" begin
        x0 = synchronized_state(N, [1.2,2.0])
        @test synchronization_error(x0) == 0
        @test lyapunov_dispersion(x0,p.tau) == 0
        @test phase_order_squared(x0,[0.0,0.0]) ≈ 1 atol=1e-14
        @test isnan(phase_order_squared(x0,[1.2,2.0]))
        @test_throws DomainError phase_order_squared(x0,[1.2,2.0]; strict=true)
        two = [1.0,0.0,-1.0,0.0]
        @test synchronization_error(two) ≈ 1.0
        @test phase_order_squared(two,[0.0,0.0]) ≈ 0.0 atol=1e-14
        @test phase_order(two,[0.0,0.0]) ≈ 0.0 atol=1e-14

        x = distributed_initial_state(N,p)
        net = FHNNetwork(N,p; lambda=-1.35)
        E = synchronization_error(x)
        V = lyapunov_dispersion(x,p.tau)
        @test N/2 * min(1,p.tau)*E <= V + 1e-10
        @test V <= N/2 * max(1,p.tau)*E + 1e-10
        @test lyapunov_derivative(x,net) <= lyapunov_upper_bound(x,net) + 1e-10
        @test_throws ArgumentError lyapunov_upper_bound(x,FHNNetwork(N,p;
                                          adjacency=ring_adjacency(N)))
        @test_throws ArgumentError lyapunov_upper_bound(x,FHNNetwork(N,p;
                                          currents=[0.1; zeros(N-1)]))
        sol = simulate_network(net; u0=x,tspan=(0.,100.),saveat=0.25)
        metr = network_metrics(sol,net)
        @test metr.E[1] ≈ E
        @test metr.V[1] ≈ V
        @test length(metr.E) == length(sol.t)
        @test minimum(metr.E) >= -1e-10
        c = min(0.35,p.b/p.tau)
        @test all(metr.V[k] <= metr.V[1]*exp(-2*c*metr.t[k]) + 1e-7
                  for k in eachindex(metr.t))
        @test all(isnan(metr.r[k]) == !metr.phase_defined[k] for k in eachindex(metr.t))
    end

    @testset "Bloques longitudinal/transversal y Floquet del ciclo" begin
        lambda = -1.5
        blocks = mean_field_linear_blocks(p,lambda)
        @test blocks.parallel ≈ jacobian_fhn(equilibrium(p)[1],p)
        @test blocks.transverse ≈ blocks.parallel + [lambda 0.0; 0.0 0.0]
        @test blocks.transverse_stable
        @test !blocks.full_equilibrium_stable  # I=1.4: equilibrio individual inestable
        @test_throws ArgumentError mean_field_linear_blocks(p,Inf)
        det = detect_limit_cycle(p; transient=300, window=300, dt=0.35)
        @test det.status == :periodic_candidate
        orb = refine_periodic_orbit(p,det)
        @test orb.status == :converged
        start = initial_state_from_orbit(orb,N)
        @test length(start) == 2N
        @test norm(start[1:2] - start[3:4]) > 1e-3
        f0 = transverse_floquet(orb,p,0.0)
        @test any(abs(z-1) < 1e-3 for z in f0.multipliers)
        fneg = transverse_floquet(orb,p,lambda)
        @test fneg.stability == :attractive
        @test fneg.spectral_radius < 1
        @test_throws ArgumentError initial_state_from_orbit(orb,1)
    end

    @testset "Gráficos, almacenamiento y separación de responsabilidades" begin
        net = FHNNetwork(N,p; lambda=-0.8)
        sol = simulate_network(net; tspan=(0.,10.), saveat=0.5)
        m = network_metrics(sol,net)
        @test size(network_state_matrix(sol,N)) == (2N,length(sol.t))
        @test length(network_layout(N)) == N
        # Guardamos datos portables, sin serializar el solucionador de EDO.
        mktempdir() do dir
            file = joinpath(dir,"prueba_C.jld2")
            @test save_network_result(file,sol,net,m) == file
            content = load_network_result(file)
            @test content["states"] ≈ network_state_matrix(sol,N)
            @test content["error_E"] ≈ m.E
            @test content["phase_r2"] ≈ m.r2
            @test content["parameters"].N == N
            @test isfile(file)
            # Prueba rápida de renderizado estático con la misma dependencia Makie.
            png = joinpath(dir,"grafo.png")
            plot_network_snapshot(sol.u[end],net; savepath=png)
            @test isfile(png) && filesize(png) > 1000
        end
    end
end
println("Etapa C: pruebas completadas.")
