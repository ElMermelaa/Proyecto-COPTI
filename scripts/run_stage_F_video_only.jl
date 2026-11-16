# Genera SOLO el MP4 comparativo de la Etapa E, sin recalcular las HJB.
# julia --project=. scripts/run_stage_F_video_only.jl
include(joinpath(@__DIR__,"..","src","FitzHughNagumoControl.jl"))
using .FitzHughNagumoControl
root=normpath(joinpath(@__DIR__,".."))
cache=joinpath(root,"data","processed","etapa_E_sync.jld2")
isfile(cache) || error("Falta etapa_E_sync.jld2; ejecuta primero scripts/run_stage_E.jl")
out=joinpath(root,"docs","videos","F_01_comparacion_control_variable_vs_constante.mp4")
println("Generando video usando trayectoria guardada de Etapa E, mismo N y x0.")
animate_stage_E_comparison(cache,out;delta=0.03,tmax=14.0,nframes=140)
println("Video generado: ",out)

cacheD=joinpath(root,"data","processed","etapa_E_desync.jld2")
if isfile(cacheD)
    videoD=joinpath(root,"docs","videos","F_02_comparacion_corrientes_vs_sin_control.mp4")
    println("Generando comparación de desincronización con el mismo N y x0.")
    animate_stage_E_desync_comparison(cacheD,videoD;nframes=140)
    println("Video generado: ",videoD)
else
    println("No existe la caché etapa_E_desync.jld2: se omite segundo video.")
end
