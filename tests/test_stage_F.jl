using Test
using LinearAlgebra
include(joinpath(@__DIR__,"..","src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

@testset "Etapa F - HJB reducida y 4D" begin
    p=FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
    c=equilibrium(p)
    @test abs(c[1]) < 1e-12

    @testset "Reducción exacta en N=2" begin
        d=[0.65,0.12];x=reduced_sync_state(d,p)
        @test length(x)==4
        @test isapprox(synchronization_error(x),sum(abs2,d);atol=1e-12)
        @test isapprox(phase_order_squared(x,c),0.0;atol=1e-12)
        for lambda in (-2.5,-1.5,0.3)
            f=sync_ocp_field(x,lambda,p)
            g=reduced_sync_rhs(d,lambda,p)
            @test isapprox(f,[g[1],g[2],-g[1],-g[2]];atol=1e-11)
        end
        @test_throws DimensionMismatch reduced_sync_state([1.0],p)
        @test_throws ArgumentError reduced_sync_state([0.5,0.1],with_current(p,0.8))
    end

    x0=reduced_sync_state([0.65,0.12],p)
    spec=SyncOCPSpec(p,2,x0;center=c,delta=0.03)
    @testset "HJB de llegada estacionaria" begin
        hs=solve_reduced_sync_hjb(spec;
             vgrid=collect(range(-0.9,0.9,length=31)),
             wgrid=collect(range(-0.45,0.45,length=25)),
             actions=[-2.5,-1.5,-0.5,0.5],
             h=0.15,maxiter=90)
        @test size(hs.W)==(31,25)
        @test hs.iterations>=1
        @test all(z>=0 for z in hs.W)
        @test reduced_sync_value(hs,[0.0,0.0])==0
        @test minimum(hs.W)==0
        @test count(isfinite,hs.W)>10
        testpoint=[0.2,0.0]
        @test isfinite(reduced_sync_value(hs,testpoint)) ||
              isinf(reduced_sync_value(hs,testpoint))
        fb=reduced_sync_feedback(hs,[0.3,0.1])
        @test spec.lambda_min<=fb.lambda<=spec.lambda_max
        rr=rollout_reduced_sync(hs,[0.25,0.05];mode=:constant,
                                lambda_constant=-1.5,dt=0.03,Tmax=15.0)
        @test rr.status==:reached
        @test rr.E[end]<=spec.delta+1e-8
        @test all(diff(rr.t).>0)
        @test issorted(rr.J)
        @test all(isfinite,rr.J)
        @test_throws ArgumentError solve_reduced_sync_hjb(
            SyncOCPSpec(p,3,distributed_initial_state(3,p);delta=1e-5))
    end

    @testset "HJB temporal completa N=2 y restricciones" begin
        sd=DesyncOCPSpec(p,2,synchronized_state(2,[1.2,1.1]);center=c,
                         T=0.4,Umax=0.8,lambda0=-0.15,
                         alpha_d=3.0,gamma_d=0.05,kappa=8.0)
        grid=solve_small_desync_hjb(sd;
             vaxis=collect(range(-2.0,2.0,length=5)),
             waxis=[0.1,0.9,1.7,2.5],actions=[-0.8,0.0,0.8],steps=2)
        @test length(grid.values)==3
        @test size(grid.values[1])==(5,4,5,4)
        @test all(z>=0 for z in grid.values[end])
        @test all(isfinite, grid.values[end])
        pol=small_desync_feedback(grid,sd.x0,1)
        @test isapprox(sum(pol.u),0.0;atol=1e-12)
        @test maximum(abs,pol.u)<=sd.Umax+1e-12
        for mode in (:hjb,:zero)
            rr=rollout_small_desync(grid;mode=mode)
            @test isfinite(rr.final_J)
            @test isfinite(rr.terminal_r2)
            @test size(rr.currents)==(2,2)
            @test maximum(abs,rr.currents)<=sd.Umax+1e-12
            @test maximum(abs.(vec(sum(rr.currents;dims=1))))<1e-12
            @test isapprox(rr.E[1],0.0;atol=1e-12)
        end
        @test_throws ArgumentError solve_small_desync_hjb(
            DesyncOCPSpec(p,3,synchronized_state(3,[1.2,1.1]);T=0.4))
    end
end
println("Etapa F: pruebas completadas.")
