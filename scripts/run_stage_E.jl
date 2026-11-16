# Desde la raíz del repositorio:
#   julia --project=. tests/test_stage_E.jl
#   julia --project=. scripts/run_stage_E.jl
# Opción avanzada de tiro indirecto (no garantiza convergencia): --shooting
# La discretización y pesos son decisiones experimentales, no datos del enunciado.

include(joinpath(@__DIR__,"..","src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl
using JLD2

root=normpath(joinpath(@__DIR__,".."))
figdir=joinpath(root,"docs","figures")
datadir=joinpath(root,"data","processed")
mkpath(figdir);mkpath(datadir)

p=FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
N=4   # Red pequeña para optimizar y diagnosticar; A-D usan también N=12.
c=equilibrium(p)
cache=joinpath(datadir,"etapa_B_resultados.jld2")

function initial_states(p,N,cache)
    if isfile(cache)
        db=JLD2.load(cache)
        q=db["parameters"]
        if all(isapprox(getproperty(p,f),getproperty(q,f);atol=1e-12)
               for f in (:a,:b,:tau,:I))
            pts=db["orbit_states"]
            x0=zeros(2N)
            for i in 1:N
                j=1+round(Int,(i-1)*(size(pts,2)-1)/N)
                x0[2i-1:2i] .= pts[:,j]
            end
            println("Inicialización: órbita refinada de Etapa B.")
            return x0,synchronized_state(N,pts[:,1])
        end
    end
    x0=distributed_initial_state(N,p)
    x0[1:2:end] .+= 0.3
    println("Inicialización alternativa: distribución determinista, no ciclo certificado.")
    return x0,synchronized_state(N,[-1.4,-0.2])
end

x0,xsync=initial_states(p,N,cache)
println("Parámetros: ",p,"; N=",N)
println("E(x0)=",synchronization_error(x0))

# SINCRONIZACIÓN: tiempo libre; caso directo de dimensión N=4.
spec=SyncOCPSpec(p,N,x0;center=c,delta=0.03,
         lambda_bounds=(-2.5,0.5),lambda_ref=-0.15,
         rho=1.0,alpha=0.1,beta=0.2,gamma=0.5)
println("\n===== SINCRONIZACIÓN =====")
rs=optimize_sync_direct(spec;M=10,substeps=8,T0=30.0,
            T_bounds=(5.0,100.0),initial_lambda=-1.35,
            outer_iterations=4,inner_iterations=35,feas_tol=5e-3)
println("T candidato = ",round(rs.T,digits=5),
        "; J_RK4 = ",round(rs.J,digits=5),
        "; E(T) discretizado = ",rs.terminal_E,
        "; E(T)-δ = ",rs.constraint_residual,
        "; factibilidad (RK4 y Tsit5) = ",rs.feasible,
        "; estado = ",rs.solver_status, "; origen = ",rs.feasible_source)
println("Control λ (por segmentos) = ",round.(rs.controls,digits=3))
println("E(T) verificado dentro del optimizador con Tsit5 = ",rs.terminal_E_independent)
rs.feasible || error("La optimización de sincronización NO consiguió E(T)<=δ. No se generarán figuras ni datos E inválidos. Prueba ampliar T_bounds o ajustar M/iteraciones.")

simS=simulate_optimal_sync(spec,rs;saveat=0.15)
actualE=synchronization_error(simS.sol.u[end])
actualE<=spec.delta || error("La segunda verificación de sincronización con Tsit5 viola E(T)<=δ. No se guardarán resultados inválidos.")
println("E(T) con Tsit5 independiente = ",actualE,
        "; discrepancia en E(T) = ",abs(actualE-rs.terminal_E))
js=evaluate_sync_cost(simS.sol,simS.network;delta=spec.delta,
       rho=spec.rho,alpha=spec.alpha,beta=spec.beta,
       gamma=spec.gamma,lambda_ref=spec.lambda_ref,center=spec.center)
println("Primera llegada a δ: ",js.status,
        "; t = ",js.arrival_time,
        "; J_llegada = ",js.arrival_J)
println("J acumulado hasta T candidato = ",js.horizon_J)

# Referencia Etapa D: λ constante -1.5 (mismo x0 y pesos).
refnet=FHNNetwork(N,p;lambda=-1.5)
refsol=simulate_network(refnet;u0=x0,tspan=(0.0,110.0),saveat=0.2)
refmet=network_metrics(refsol,refnet;center=c)
refcost=evaluate_sync_cost(refsol,refnet;delta=spec.delta,
     rho=spec.rho,alpha=spec.alpha,beta=spec.beta,
     gamma=spec.gamma,lambda_ref=spec.lambda_ref,center=c)
println("Comparación λ=-1.5 (Etapa D): llegada=",refcost.status,
        ", T=",refcost.arrival_time,", J_llegada=",refcost.arrival_J)

diagS=pmp_diagnostics_sync(spec,rs,simS.sol;samples=120)
println("Residuo máximo de la ley λ_PMP: ",diagS.max_projected_control_residual)
println("Residuo H_s(T)=0: ",diagS.terminal_H,
        "; μ candidato: ",diagS.mu)
plot_optimal_sync(rs,spec,simS,js;
     comparison=(label="λ=-1.5",metrics=refmet),
     savepath=joinpath(figdir,"E_01_sincronizacion_candidata.png"))
plot_pmp_sync(diagS;savepath=joinpath(figdir,"E_02_residuos_PMP_sync.png"))
save_optimal_result(joinpath(datadir,"etapa_E_sync.jld2"),spec,rs,simS,diagS)

# DESINCRONIZACIÓN: corriente con suma cero, T fijo, multistart.
println("\n===== DESINCRONIZACIÓN =====")
specD=DesyncOCPSpec(p,N,xsync;center=c,T=28.0,lambda0=-0.15,
             Umax=0.8,alpha_d=3.0,gamma_d=0.05,kappa=8.0)
rd=optimize_desync_direct(specD;M=10,substeps=8,
              maxiter=30,seed_amplitude=0.65)
println("J candidato = ",rd.J,"; J sin corriente = ",rd.baseline_J)
println("r²(T) discretizado = ",rd.terminal_r2,
        "; max |u_i| = ",rd.max_current,
        "; max |Σu_i| = ",rd.max_current_sum,
        "; semilla ganadora = ",rd.best_seed)

simD=simulate_optimal_desync(specD,rd;saveat=0.15)
jd=evaluate_desync_cost(simD.sol,simD.network;
     alpha_d=specD.alpha_d,gamma_d=specD.gamma_d,
     kappa=specD.kappa,Umax=specD.Umax,center=specD.center)
println("J verificado Tsit5 = ",jd.total_J,
        "; r²(T) = ",jd.terminal_r2,
        "; E(T) = ",jd.terminal_E)
diagD=pmp_diagnostics_desync(specD,rd,simD.sol;samples=120)
println("Residuo máximo de ley u_PMP: ",diagD.max_projected_control_residual)
plot_optimal_desync(rd,specD,simD,jd;
     savepath=joinpath(figdir,"E_03_desincronizacion_candidata.png"))
plot_pmp_desync(diagD;savepath=joinpath(figdir,"E_04_residuos_PMP_desync.png"))
save_optimal_result(joinpath(datadir,"etapa_E_desync.jld2"),specD,rd,simD,diagD)

if "--shooting" in ARGS
    println("\n===== SHOOTING INDIRECTO (experimental) =====")
    shoot=shoot_sync_pmp(spec;p0=diagS.adjoint_solution(0.0),
                           T0=rs.T,mu0=max(rs.mu,1e-3),maxiters=20)
    println(shoot)
end
println("\nFiguras guardadas en: ",figdir)
println("Resultados guardados en: ",datadir)
println("Etapa E: revisar factibilidad y residuos antes de hablar de optimalidad.")
