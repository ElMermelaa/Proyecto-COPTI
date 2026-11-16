# Persistencia portable de candidatos y métricas. Sin ODESolution ni closures.

function save_optimal_result(path,s::Union{SyncOCPSpec,DesyncOCPSpec},result,sim,diag;
                             description="Etapa E: candidato de control óptimo")
    file=abspath(string(path));mkpath(dirname(file))
    p=s.neuron
    pars=(a=p.a,b=p.b,tau=p.tau,I=p.I,N=s.N)
    T=Float64(result.T)
    controls=Float64.(result.controls)
    λ_samples=Float64[_law_at(sim.network.lambda,t) for t in sim.sol.t]
    U_samples=Float64.(control_current_matrix(sim.sol,sim.network))
    JLD2.jldsave(file;
        format_version="FHN-stage-E-1",saved_at=string(Dates.now()),
        description=string(description),kind=string(result.kind),
        parameters=pars,initial_state=Float64.(s.x0),center=Float64.(s.center),
        T=T,control_segments=controls,
        segment_boundaries=collect(range(0.0,T,length=result.M+1)),
        times=Float64.(sim.sol.t),states=network_state_matrix(sim.sol,s.N),
        error_E=sim.metrics.E,phase_r2=sim.metrics.r2,
        lambda_samples=λ_samples,current_samples=U_samples,
        objective_RK4=Float64(result.J),
        terminal_E=Float64(synchronization_error(sim.sol.u[end])),
        pmp_control_residual=Float64(diag.max_projected_control_residual),
        pmp_terminal_residual=Float64(diag.terminal_adjoint_residual),
        algorithm=string(result.method),M=result.M,substeps=result.substeps)
    return file
end

function load_optimal_result(path)
    data=JLD2.load(string(path))
    get(data,"format_version",nothing)=="FHN-stage-E-1" ||
        throw(ArgumentError("Archivo incompatible con Etapa E."))
    return data
end
