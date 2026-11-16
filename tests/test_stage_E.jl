using Test
using LinearAlgebra
using ForwardDiff
using JLD2

include(joinpath(@__DIR__,"..","src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa E - Pontryagin y controles candidatos" begin
    p=FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
    N=4
    c=equilibrium(p)
    x=distributed_initial_state(N,p;radius_v=1.2,radius_w=0.3)
    x[1:2:end] .+= 0.35 # evitar promedio en equilibrio inestable
    ss=SyncOCPSpec(p,N,x;delta=0.15,lambda_bounds=(-2.0,0.0),
                   lambda_ref=-0.15,rho=1.0,alpha=0.1,beta=0.2,gamma=0.5)
    zs=synchronized_state(N,[-1.2,-0.1])
    sd=DesyncOCPSpec(p,N,zs;T=4.0,lambda0=-0.15,Umax=0.8,
                     alpha_d=1.5,gamma_d=0.2,kappa=3.0)

    @testset "Datos de los problemas" begin
        @test ss.N==N && sd.N==N
        @test_throws ArgumentError SyncOCPSpec(p,N,zs;delta=0.01)
        @test_throws ArgumentError SyncOCPSpec(p,N,x;delta=0.0)
        @test_throws ArgumentError SyncOCPSpec(p,N,x;lambda_bounds=(0.0,0.0))
        @test_throws ArgumentError DesyncOCPSpec(p,N,zs;Umax=0.0)
        @test_throws DimensionMismatch DesyncOCPSpec(p,N,zeros(3))
    end

    @testset "Gradientes y Hamiltonianos por diferenciación automática" begin
        gE=sync_error_gradient(x)
        gR=phase_order_squared_gradient(x,c)
        @test gE ≈ ForwardDiff.gradient(y->FitzHughNagumoControl._ocp_E(y),x) atol=1e-10
        @test gR ≈ ForwardDiff.gradient(y->FitzHughNagumoControl._ocp_r2(y,c),x) atol=1e-8
        @test abs(sum(gE[1:2:end]))<1e-12
        @test abs(sum(gE[2:2:end]))<1e-12
        @test_throws DomainError phase_order_squared_gradient(synchronized_state(N,c),c)
        λ=-0.4
        @test sync_ocp_field(x,λ,p) ≈ network_rhs(x,FHNNetwork(N,p;lambda=λ)) atol=1e-12
        pv=collect(range(-0.4,0.6,length=2N))
        @test sync_adjoint(x,pv,λ,ss) ≈
              -ForwardDiff.gradient(y->sync_hamiltonian(y,pv,λ,ss),x) atol=2e-8
        u=[0.2,-0.2,0.1,-0.1]
        @test desync_ocp_field(x,u,-0.15,p) ≈
              network_rhs(x,FHNNetwork(N,p;lambda=-0.15,currents=u)) atol=1e-12
        @test desync_adjoint(x,pv,u,sd) ≈
              -ForwardDiff.gradient(y->desync_hamiltonian(y,pv,u,sd),x) atol=2e-8
    end

    @testset "Mínimos puntuales y proyecciones" begin
        pv=collect(range(-2.0,1.5,length=2N))
        λstar=sync_control_pmp(x,pv,ss)
        @test ss.lambda_min<=λstar<=ss.lambda_max
        @test all(sync_hamiltonian(x,pv,λstar,ss) <=
                  sync_hamiltonian(x,pv,λ,ss)+1e-10
                  for λ in range(ss.lambda_min,ss.lambda_max,length=31))
        Ustar=desync_control_pmp(pv,sd)
        @test abs(sum(Ustar))<1e-9
        @test maximum(abs,Ustar)<=sd.Umax+1e-12
        @test desync_hamiltonian(x,pv,Ustar,sd)<=
              desync_hamiltonian(x,pv,zeros(N),sd)+1e-10
        @test project_balanced_box([2.0,-1.0,-1.0,0.0],0.5) ≈
              [0.5,-0.5,-0.5,0.5] atol=1e-9
        @test_throws ArgumentError project_balanced_box(ones(3),-1)
    end

    @testset "Discretización y optimizadores directos" begin
        @test piecewise_value([-1.0,-0.5],0.0,4.0)==-1.0
        @test piecewise_value([-1.0,-0.5],4.0,4.0)==-0.5
        @test piecewise_currents([0.1 0.3; -0.1 -0.3],3.9,4.0) ≈ [0.3,-0.3]
        result=optimize_sync_direct(ss;M=5,substeps=4,T0=8.0,
             T_bounds=(3.0,12.0),initial_lambda=-1.3,
             penalty0=3.0,outer_iterations=2,inner_iterations=5)
        @test length(result.controls)==5
        @test 3.0<=result.T<=12.0
        @test all(ss.lambda_min-1e-12<=v<=ss.lambda_max+1e-12 for v in result.controls)
        @test result.J>=ss.rho*result.T-1e-7
        @test isfinite(result.terminal_E)
        @test abs(result.constraint_residual-(result.terminal_E-ss.delta))<1e-10
        @test result.feasible == (result.constraint_residual<=0.0 &&
                isfinite(result.terminal_E_independent) &&
                result.terminal_E_independent<=ss.delta)

        rd=optimize_desync_direct(sd;M=4,substeps=4,maxiter=4,seed_amplitude=0.4)
        @test size(rd.controls)==(N,4)
        @test all(abs(sum(rd.controls[:,j]))<1e-9 for j in 1:4)
        @test maximum(abs,rd.controls)<=sd.Umax+1e-10
        @test rd.J<=rd.baseline_J+1e-8
        @test length(rd.seed_costs)>=3
        @test isfinite(rd.terminal_r2)
        @test 0<=rd.terminal_r2<=1+1e-10
    end

    @testset "Soluciones de EDO, adjuntos y persistencia" begin
        # Primero validamos candidatos prefijados, sin asumir que sean óptimos.
        rS=(controls=fill(-1.4,4),M=4,T=3.0,mu=1.0)
        rD=(controls=zeros(N,4),M=4,T=sd.T)
        simS=simulate_optimal_sync(ss,rS;saveat=0.25)
        simD=simulate_optimal_desync(sd,rD;saveat=0.25)
        diagS=pmp_diagnostics_sync(ss,rS,simS.sol;samples=12)
        diagD=pmp_diagnostics_desync(sd,rD,simD.sol;samples=12)
        @test length(simS.sol.u[end])==2N
        @test length(diagS.lambda_pmp)==12
        @test isfinite(diagS.terminal_H)
        @test diagS.terminal_adjoint_residual<1e-7
        @test diagD.terminal_adjoint_residual<1e-7
        @test diagD.max_current_sum<1e-12
        mktempdir() do dir
            file=save_optimal_result(joinpath(dir,"optimal.jld2"),ss,
                     merge(rS,(kind=:synchronization,J=1.0,substeps=4,
                     method=:testing)),simS,diagS)
            @test isfile(file)
            data=load_optimal_result(file)
            @test data["format_version"]=="FHN-stage-E-1"
            @test size(data["states"],1)==2N
            plot_pmp_sync(diagS;savepath=joinpath(dir,"pmp.png"))
            @test isfile(joinpath(dir,"pmp.png"))
        end
    end
    @testset "Restricción dura y verificación Tsit5" begin
        xnear = synchronized_state(2, [1.0, 0.7])
        xnear[1] += 0.3
        sc = SyncOCPSpec(p,2,xnear; delta=0.005,
                lambda_bounds=(-2.5,0.0),lambda_ref=-0.15)
        rc=optimize_sync_direct(sc;M=3,substeps=6,T0=6.0,
                T_bounds=(1.0,90.0),initial_lambda=-1.5,
                outer_iterations=1,inner_iterations=2)
        @test rc.feasible
        @test rc.terminal_E<=sc.delta
        @test rc.terminal_E_independent<=sc.delta
        simc=simulate_optimal_sync(sc,rc;saveat=0.2)
        @test synchronization_error(simc.sol.u[end])<=sc.delta
        @test all(sc.lambda_min<=z<=sc.lambda_max for z in rc.controls)
    end

end
println("Etapa E: tests completados.")
