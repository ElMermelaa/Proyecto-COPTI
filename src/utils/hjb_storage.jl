# Guardado de matrices y metadatos SIN closures ni soluciones de EDO.
function save_hjb_result(file,ss::ReducedSyncHJB,sd::SmallDesyncHJB,
        sh,sc,dh,dz)
    path=abspath(string(file));mkpath(dirname(path))
    JLD2.jldsave(path;
        format_version="FHN-stage-F-1",
        description="HJB exploratoria: sincronización reducida exacta y desincronización N=2 4D",
        sync_vgrid=ss.vgrid,sync_wgrid=ss.wgrid,sync_W=ss.W,
        sync_policy=ss.policy,sync_actions=ss.actions,sync_delta=ss.spec.delta,
        sync_iterations=ss.iterations,
        sync_hjb_time=sh.t,sync_hjb_E=sh.E,sync_hjb_J=sh.J[end],
        sync_hjb_lambda=sh.lambda,sync_hjb_uncovered=sh.uncovered_steps,
        sync_constant_time=sc.t,sync_constant_E=sc.E,sync_constant_J=sc.J[end],
        desync_axes=sd.axes,desync_actions=sd.actions,desync_timestep=sd.h,
        desync_value_time0=sd.values[1],
        desync_hjb_t=dh.t,desync_hjb_r2=dh.r2,desync_hjb_J=dh.final_J,
        desync_hjb_covered=dh.covered,desync_hjb_currents=dh.currents,
        desync_zero_r2=dz.r2,desync_zero_J=dz.final_J)
    return path
end
