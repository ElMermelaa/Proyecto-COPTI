# ETAPA D — Datos portables, sin serializar funciones ni ODESolution.

"""    save_control_result(path, sol, system, metrics;
                               label="", kind="", sync_cost=nothing,
                               desync_cost=nothing, reltol=1e-9, abstol=1e-9)

Guarda parámetros, estados, métricas, acoplamientos y controles efectivamente
aplicados. Los cierres no se serializan; las corrientes u_i(t) quedan muestreadas.
"""
function save_control_result(path,sol,
        system::Union{FHNNetwork,FHNCurrentFeedback},metrics;
        label="",kind="",sync_cost=nothing,desync_cost=nothing,
        reltol::Real=1e-9,abstol::Real=1e-9,description="Etapa D")
    net = system isa FHNNetwork ? system : system.base
    file = abspath(string(path))
    mkpath(dirname(file))
    p = net.neuron
    params = (a=p.a,b=p.b,tau=p.tau,I=p.I,N=net.N)
    gain = system isa FHNCurrentFeedback ? system.k : 0.0
    lambdas = Float64[_law_at(net.lambda,t) for t in sol.t]
    U = control_current_matrix(sol,system)
    j_sync = sync_cost === nothing ? NaN : sync_cost.arrival_J
    j_desync = desync_cost === nothing ? NaN : desync_cost.total_J
    arrival_time = sync_cost === nothing ? NaN : sync_cost.arrival_time
    status = sync_cost === nothing ? "not_evaluated" : string(sync_cost.status)
    JLD2.jldsave(file;
        format_version="FHN-stage-D-1",
        description=string(description), saved_at=string(Dates.now()),
        label=string(label),kind=string(kind),
        parameters=params,adjacency=net.A,
        initial_state=Float64.(sol.u[1]),
        times=Float64.(sol.t),states=network_state_matrix(sol,net.N),
        center=metrics.center,
        error_E=metrics.E,lyapunov_V=metrics.V,
        phase_r=metrics.r,phase_r2=metrics.r2,
        min_center_distance=metrics.min_center_distance,
        phase_defined=metrics.phase_defined,
        lambda_samples=lambdas,current_samples=U,
        feedback_gain=gain,
        sync_arrival_J=j_sync,sync_arrival_time=arrival_time,
        sync_arrival_status=status,desync_total_J=j_desync,
        integrator="Tsit5",reltol=Float64(reltol),abstol=Float64(abstol))
    return file
end

"""    load_control_result(path)

Lee resultados planos de la Etapa D para analizar o graficar sin nueva EDO.
"""
function load_control_result(path)
    data=JLD2.load(string(path))
    get(data,"format_version",nothing) == "FHN-stage-D-1" ||
        throw(ArgumentError("El archivo no pertenece a la Etapa D."))
    return data
end
