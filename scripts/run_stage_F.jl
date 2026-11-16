# Desde la raíz del repo:
# julia --project=. tests/test_stage_F.jl
# julia --project=. scripts/run_stage_F.jl
# julia --project=. scripts/run_stage_F.jl --video
# El video compara DOS controles de Etapa E para N=4; NO compara HJB y PMP.

include(joinpath(@__DIR__,"..","src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl

root=normpath(joinpath(@__DIR__,".."))
figdir=joinpath(root,"docs","figures")
viddir=joinpath(root,"docs","videos")
datadir=joinpath(root,"data","processed")
mkpath(figdir);mkpath(viddir);mkpath(datadir)

p=FHNParameters(a=0.7,b=0.5,tau=12.5,I=1.4)
center=equilibrium(p)
println("Etapa F: equilibrio c=",center)
println("NOTA: reducción de sincronización exacta SOLO en N=2, v*=0.")

# HJB estacionaria de sincronización en subvariedad exacta N=2.
d0=[0.65,0.12]
x0=reduced_sync_state(d0,p)
spec=SyncOCPSpec(p,2,x0;center=center,delta=0.03,
            lambda_bounds=(-2.5,0.5),lambda_ref=-0.15,
            rho=1.0,alpha=0.1,beta=0.2,gamma=0.5)
println("\n===== HJB SINCRONIZACIÓN REDUCIDA =====")
hs=solve_reduced_sync_hjb(spec)
valid=count(isfinite,hs.W)
println("Celdas con valor finito: $valid / $(length(hs.W))")
println("Iteraciones Bellman: ",hs.iterations,
        "; valor aproximado en d0: ",reduced_sync_value(hs,d0))
sh=rollout_reduced_sync(hs,d0;mode=:hjb,dt=0.05,Tmax=25.0)
sc=rollout_reduced_sync(hs,d0;mode=:constant,lambda_constant=-1.5,dt=0.05,Tmax=25.0)
println("Feedback HJB: estado=",sh.status,", T=",sh.t[end],", J=",sh.J[end],
        ", pasos fuera de cobertura=",sh.uncovered_steps)
println("Referencia λ=-1.5: estado=",sc.status,", T=",sc.t[end],", J=",sc.J[end])
println("La reducción obliga a r²=0 fuera del objetivo; no representa sincronización de fase general.")
plot_hjb_sync(hs;d0=d0,savepath=joinpath(figdir,"F_01_hjb_sync_reducida.png"))

# HJB temporal para desincronización N=2 COMPLETO (4 dimensiones).
println("\n===== HJB DESINCRONIZACIÓN N=2 =====")
xsync=synchronized_state(2,[1.2,1.1])
specD=DesyncOCPSpec(p,2,xsync;center=center,T=1.2,lambda0=-0.15,
                    Umax=0.8,alpha_d=3.0,gamma_d=0.05,kappa=8.0)
hd=solve_small_desync_hjb(specD;steps=6)
dh=rollout_small_desync(hd;mode=:hjb)
dz=rollout_small_desync(hd;mode=:zero)
println("HJB discreta: J=",dh.final_J,", r²(T)=",dh.terminal_r2,
        ", dentro de grilla=",dh.covered,", control máximo=",maximum(abs,dh.currents))
if !dh.covered
    println("ADVERTENCIA: fuera de cobertura HJB. Se empleo control u=0 como " *
            "fallback en celdas sin valores finitos. No interpretar como " *
            "politica HJB valida ni como optimo aproximado en toda la trayectoria.")
end

println("Referencia u=0: J=",dz.final_J,", r²(T)=",dz.terminal_r2)
println("Horizonte 1.2, malla espacial gruesa y acciones discretas: resultados EXPLORATORIOS.")
plot_hjb_desync(hd;savepath=joinpath(figdir,"F_02_hjb_desync_corte_4D.png"))
plot_hjb_comparison(sh,sc,dh,dz;savepath=joinpath(figdir,"F_03_comparaciones_HJB.png"))

out=save_hjb_result(joinpath(datadir,"etapa_F_hjb.jld2"),hs,hd,sh,sc,dh,dz)
println("Figuras: ",figdir)
println("Datos: ",out)

if "--video" in ARGS
    cache=joinpath(datadir,"etapa_E_sync.jld2")
    if isfile(cache)
        path=joinpath(viddir,"F_01_comparacion_control_variable_vs_constante.mp4")
        animate_stage_E_comparison(cache,path;delta=0.03,tmax=14.0,nframes=140)
        println("Video comparativo Etapa E generado: ",path)
        cacheD=joinpath(datadir,"etapa_E_desync.jld2")
        if isfile(cacheD)
            pathD=joinpath(viddir,"F_02_comparacion_corrientes_vs_sin_control.mp4")
            animate_stage_E_desync_comparison(cacheD,pathD;nframes=140)
            println("Video comparativo de desincronización: ",pathD)
        end
    else
        println("No se encuentra etapa_E_sync.jld2: no se genera video; ejecuta primero run_stage_E.jl.")
    end
end
println("Etapa F terminada. No implica solución global certificada de HJB.")
