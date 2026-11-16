using Test
using LinearAlgebra
using JLD2

include(joinpath(@__DIR__, "..", "src", "FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa D — controles suficientes y costos" begin
    p = FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
    N = 6
    A = complete_adjacency(N)
    x = distributed_initial_state(N,p)
    center = equilibrium(p)

    @testset "Realimentación: estructura y equivalencia de campo medio" begin
        U = feedback_currents(x,1.4)
        @test length(U) == N
        @test abs(sum(U)) < 1e-12
        @test U ≈ [-1.4*(x[2i-1]-sum(x[1:2:end])/N) for i in 1:N]
        @test_throws ArgumentError feedback_currents(x,-0.1)
        @test_throws ArgumentError feedback_currents(x,Inf)
        net = FHNNetwork(N,p;adjacency=A,lambda=-0.2)
        ctrl = FHNCurrentFeedback(net;k=1.4)
        rhs = feedback_rhs(x,ctrl)
        # k y lambda son algebraicamente intercambiables SOLO en campo medio.
        equivalent = FHNNetwork(N,p;adjacency=A,lambda=-1.6)
        @test rhs ≈ network_rhs(x,equivalent) atol=1e-12
        inplace = similar(x)
        feedback_rhs!(inplace,x,ctrl,0.0)
        @test inplace ≈ rhs atol=1e-12
        @test_throws ArgumentError FHNCurrentFeedback(net;k=-2)
        @test_throws DimensionMismatch feedback_rhs(zeros(4),ctrl)
        @test feedback_lyapunov_derivative(x,ctrl) <=
              feedback_lyapunov_upper_bound(x,ctrl)+1e-10
        @test_throws ArgumentError feedback_lyapunov_upper_bound(x,
                FHNCurrentFeedback(FHNNetwork(N,p;
                  adjacency=ring_adjacency(N),lambda=-0.2);k=1.4))
    end

    @testset "Criterios suficientes y tiempo garantizado" begin
        c1 = sufficient_decay_rate(p;lambda=-1.5)
        @test c1 ≈ min(0.5,p.b/p.tau)
        c2 = sufficient_decay_rate(p;lambda=0.0,k=1.5)
        @test c2 ≈ c1
        @test_throws ArgumentError sufficient_decay_rate(p;lambda=-1.0)
        @test_throws ArgumentError sufficient_decay_rate(p;lambda=-0.15)
        @test_throws ArgumentError sufficient_decay_rate(p;lambda=0,k=1)
        @test_throws ArgumentError sufficient_decay_rate(p;lambda=0,k=-1)
        p_b0 = FHNParameters(a=0.7,b=0.0,tau=12.5,I=1.4)
        @test_throws ArgumentError sufficient_decay_rate(p_b0;lambda=-2)
        E0 = synchronization_error(x)
        delta = E0/10
        Ctau = max(1,p.tau)/min(1,p.tau)
        expected = max(0.0,log(Ctau*E0/delta)/(2c1))
        @test sufficient_arrival_time(E0,delta,p;lambda=-1.5) ≈ expected
        @test sufficient_arrival_time(0.0,delta,p;lambda=-1.5) == 0.0
        @test_throws ArgumentError sufficient_arrival_time(E0,0,p;lambda=-1.5)
    end

    @testset "Perfiles balanceados y ruptura inicial de sincronía" begin
        for M in (5,6,7,12), typ in (:split,:alternating,:sinusoidal)
            pat = balanced_current_profile(M;pattern=typ)
            @test length(pat) == M
            @test abs(sum(pat)) < 1e-12
            @test isapprox(maximum(abs,pat),1.0;atol=1e-12)
            law = balanced_current_law(pat;amplitude=0.4)
            @test abs(sum(law(0.0))) < 1e-12
            @test maximum(abs,law(4.0)) <= 0.4000000001
            @test abs(sum(law(4.0))) < 1e-12
        end
        @test_throws ArgumentError balanced_current_profile(1)
        @test_throws ArgumentError balanced_current_profile(6;pattern=:unknown)
        @test_throws ArgumentError balanced_current_law([1.0,0.0])
        @test_throws ArgumentError balanced_current_law([1.0,-1.0];amplitude=-1)
        @test_throws ArgumentError balanced_current_law([1.0,-1.0];ramp_time=-1)
        smooth = balanced_current_law([1.0,-1.0];amplitude=0.8,ramp_time=1.5)
        @test all(smooth(0.0) .== 0)
        @test maximum(abs,smooth(3.0)) < 0.8

        z = [-1.1,0.0]
        x_sync = synchronized_state(N,z)
        pulse = 0.3 .* balanced_current_profile(N;pattern=:split)
        forced = FHNNetwork(N,p;lambda=-0.15,currents=pulse)
        base = FHNNetwork(N,p;lambda=-0.15)
        df = network_rhs(x_sync,forced) - network_rhs(x_sync,base)
        @test df[1:2:end] ≈ pulse
        @test all(df[2:2:end] .== 0)
        sol = simulate_network(forced;u0=x_sync,tspan=(0.0,8.0),saveat=0.1)
        @test synchronization_error(sol.u[1]) < 1e-25
        @test synchronization_error(sol.u[end]) > 1e-6
        m = network_metrics(sol,forced;center=center)
        jd = evaluate_desync_cost(sol,forced;Umax=0.5,center=center)
        @test jd.total_J > jd.running_J > 0
        @test jd.max_current <= 0.5
        @test jd.max_current_sum < 1e-10
        @test jd.initial_E < 1e-25
        @test jd.terminal_E > 0
        @test jd.terminal_r2 >= 0 && jd.terminal_r2 <= 1
        @test jd.total_J ≈ jd.phase_running_J+jd.current_running_J+jd.terminal_J atol=1e-7
        @test_throws ArgumentError evaluate_desync_cost(sol,forced;Umax=0.1)
        @test_throws ArgumentError evaluate_desync_cost(sol,
              FHNNetwork(N,p;lambda=-0.15,currents=[0.1;zeros(N-1)]))
        @test_throws ArgumentError evaluate_desync_cost(sol,
              FHNNetwork(N,p;lambda=t->-0.15))
        @test_throws ArgumentError evaluate_sync_cost(sol,forced;center=center)
    end

    @testset "Simulación controlada, Lyapunov y costos J_s" begin
        net = FHNNetwork(N,p;lambda=-1.5)
        fb = FHNCurrentFeedback(FHNNetwork(N,p;lambda=0.0);k=1.5)
        # distributed_initial_state tiene promedio casi exactamente igual al
        # equilibrio INESTABLE (v*,w*)=(0,1.4). Errores de redondeo minúsculos
        # en ese modo común se amplifican durante t=40 y pueden separar las
        # soluciones aunque ambos campos sean algebraicamente idénticos.
        # Desplazamos todas las v_i por igual: E(0) se conserva y el ensayo
        # compara la misma dinámica lejos de ese equilibrio sensible.
        x_compare = copy(x)
        x_compare[1:2:end] .+= 0.4
        sol = simulate_network(net;u0=x_compare,tspan=(0.0,40.0),saveat=0.2)
        solfb = simulate_feedback(fb;u0=x_compare,tspan=(0.0,40.0),saveat=0.2)
        @test length(sol.t)==length(solfb.t)
        @test maximum(norm(sol.u[j]-solfb.u[j]) for j in eachindex(sol.t)) < 1e-6
        cert = lyapunov_certificate(sol,net)
        certfb = lyapunov_certificate(solfb,fb)
        @test cert.c ≈ min(0.5,p.b/p.tau)
        @test certfb.c ≈ cert.c
        @test cert.max_V_violation < 1e-6
        @test cert.max_E_violation < 1e-6
        @test certfb.max_V_violation < 1e-6
        @test_throws ArgumentError lyapunov_certificate(sol,
              FHNNetwork(N,p;lambda=-0.15))
        @test_throws ArgumentError lyapunov_certificate(sol,
              FHNNetwork(N,p;lambda=t->-1.5))
        modulated = AttractiveModulatedCoupling(-1.35,0.3,0.12)
        @test all(modulated(t) <= -1.35 for t in (0.0,1.0,10.0,100.0))
        modnet = FHNNetwork(N,p;lambda=modulated)
        modsol = simulate_network(modnet;u0=x,tspan=(0.0,35.0),saveat=0.25)
        modcert = lyapunov_certificate(modsol,modnet)
        @test modcert.certified_lambda_upper == -1.35
        @test modcert.c ≈ min(0.35,p.b/p.tau)
        @test modcert.max_V_violation < 1e-6
        @test_throws ArgumentError AttractiveModulatedCoupling(-1.5,-0.1,0.2)
        @test_throws ArgumentError lyapunov_certificate(sol,
              FHNNetwork(N,p;adjacency=ring_adjacency(N),lambda=-1.5))
        @test_throws ArgumentError lyapunov_certificate(sol,
              FHNNetwork(N,p;lambda=-1.5,currents=t->zeros(N)))

        @test size(control_current_matrix(solfb,fb))==(N,length(solfb.t))
        @test total_current_effort(solfb,fb)>0
        @test total_current_effort(sol,net)≈0.0
        cost = evaluate_sync_cost(sol,net;delta=0.01,center=center)
        @test cost.horizon_J > 0
        @test cost.arrival_J >= 0
        @test cost.status in (:reached,:already_reached,:not_reached)
        @test cost.horizon_J ≈ sum(values(cost.horizon_components)) atol=1e-6
        @test first(cost.accumulated_cost) == 0
        @test last(cost.accumulated_cost) ≈ cost.horizon_J
        @test isinf(cost.arrival_J) == (!cost.reached)
        hit = first_tolerance_hit(sol,0.01)
        @test hit.status == cost.status
        if cost.reached
            @test abs(hit.E-0.01)<1e-8
            @test cost.arrival_time <= last(sol.t)
        end
        @test_throws ArgumentError evaluate_sync_cost(sol,net;rho=-1)
        @test_throws ArgumentError first_tolerance_hit(sol,0)
    end

    @testset "Casos en sincronía exacta, fase indefinida y archivos" begin
        net = FHNNetwork(N,p;lambda=-1.5)
        z = [-1.4,-0.2]
        synced = synchronized_state(N,z)
        sol = simulate_network(net;u0=synced,tspan=(0.0,5.0),saveat=0.1)
        m = network_metrics(sol,net;center=center)
        delta=0.01
        hit = first_tolerance_hit(sol,delta)
        @test hit.status == :already_reached
        @test hit.time == 0.0
        cs = evaluate_sync_cost(sol,net;delta=delta,rho=1.0,alpha=0.1,
                                beta=0.2,gamma=0.5,lambda_ref=-1.5,center=center)
        @test cs.arrival_J==0.0
        @test cs.arrival_time==0.0
        @test abs(cs.horizon_J-5.0) < 1e-5
        # Las fases no existen cuando las neuronas coinciden con el centro.
        undef = simulate_network(net;u0=synchronized_state(N,center),
                                 tspan=(0.0,1.0),saveat=0.1)
        @test_throws DomainError evaluate_sync_cost(undef,net;center=center)
        @test_throws DomainError evaluate_desync_cost(undef,net;center=center)
        @testset "Persistencia y gráficos" begin
            mktempdir() do dir
                cache = joinpath(dir,"casoD.jld2")
                @test save_control_result(cache,sol,net,m;
                           sync_cost=cs,label="prueba",kind="coupling") == cache
                loaded = load_control_result(cache)
                @test loaded["format_version"] == "FHN-stage-D-1"
                @test loaded["states"] ≈ network_state_matrix(sol,N)
                @test loaded["sync_arrival_J"] == 0.0
                @test size(loaded["current_samples"]) == (N,length(sol.t))
                @test loaded["parameters"].N == N
                cert = lyapunov_certificate(sol,net)
                file1 = joinpath(dir,"cota.png")
                file2 = joinpath(dir,"controles.png")
                plot_lyapunov_certificate(cert;savepath=file1)
                plot_control_inputs(sol,net;savepath=file2)
                @test isfile(file1) && filesize(file1)>1000
                @test isfile(file2) && filesize(file2)>1000
            end
        end
    end
end
println("Etapa D: pruebas completadas.")
