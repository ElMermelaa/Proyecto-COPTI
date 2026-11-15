# Persistencia portátil. No serializar cierres/ODEProblem/ODESolution en JLD2.

"""    save_network_result(path, sol, net, metrics; solver="Tsit5", ...)

Guarda matrices y metadatos legibles por JLD2. Permite regenerar figuras sin
repetir la simulación. Las funciones de control variables se describen pero
NO se serializan: para reproducibilidad deben documentarse en el script.
"""
function save_network_result(path, sol, net::FHNNetwork, metrics;
        solver="Tsit5", reltol=1e-9, abstol=1e-9,
        description="Simulación Etapa C")
    filename = abspath(string(path))
    mkpath(dirname(filename))
    pars = (a=net.neuron.a, b=net.neuron.b, tau=net.neuron.tau,
            I=net.neuron.I, N=net.N)
    lambda_law = net.lambda isa Real ? "constante $(net.lambda)" : "lambda(t): función del script"
    current_law = net.currents isa AbstractVector ?
        "constante individual" : "corrientes(t): función del script"
    current_vector = net.currents isa AbstractVector ? copy(net.currents) : Float64[]
    JLD2.jldsave(filename;
        format_version="FHN-stage-C-1",
        description=string(description),
        saved_at=string(Dates.now()),
        parameters=pars,
        adjacency=net.A,
        initial_state=Float64.(sol.u[1]),
        times=Float64.(sol.t),
        states=network_state_matrix(sol, net.N),
        center=metrics.center,
        error_E=metrics.E,
        lyapunov_V=metrics.V,
        phase_r=metrics.r,
        phase_r2=metrics.r2,
        phase_defined=metrics.phase_defined,
        center_distances=metrics.min_center_distance,
        lambda_description=lambda_law,
        currents_description=current_law,
        currents_constant=current_vector,
        integrator=string(solver), reltol=Float64(reltol), abstol=Float64(abstol))
    return filename
end

"""    load_network_result(path)

Devuelve `Dict` de datos planos. No requiere la solución SciML original.
"""
function load_network_result(path)
    data = JLD2.load(string(path))
    get(data, "format_version", nothing) == "FHN-stage-C-1" ||
        throw(ArgumentError("El archivo no es un resultado de la Etapa C."))
    return data
end
