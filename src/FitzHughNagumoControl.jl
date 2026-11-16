module FitzHughNagumoControl

using LinearAlgebra
using OrdinaryDiffEq
using CairoMakie
using JLD2
using Dates
using ForwardDiff
import NonlinearSolve

# Etapa A — se conservan todas las firmas anteriores.
include("models/parameters.jl")
include("models/individual.jl")
include("analysis/equilibria.jl")
include("analysis/stability.jl")
include("visualization/plots.jl")

# Etapa B — sin cambios funcionales.
include("analysis/bifurcations.jl")
include("visualization/bifurcation_plots.jl")

# Etapa C — redes, métricas, estabilidad transversal y visualización.
include("models/network.jl")
include("analysis/network_metrics.jl")
include("analysis/network_stability.jl")
include("visualization/network_plots.jl")
include("visualization/network_animations.jl")
include("utils/network_storage.jl")

export FHNParameters, with_current, with_parameter
export fhn_nonlinearity, fhn_vector_field, fhn_rhs, fhn!, nullcline_v, nullcline_w
export simulate_individual, simulate_ensemble
export equilibrium, equilibrium_residual
export jacobian_fhn, stability_info, equilibrium_sweep, parameter_sweep, hopf_currents
export plot_timeseries, plot_phase_portrait, plot_nullclines, plot_equilibrium_sweep

export hopf_analysis, upward_crossings, detect_limit_cycle
export refine_periodic_orbit, continue_periodic_branch, hopf_seed
export plot_hopf_diagram, plot_periodic_orbit, plot_cycle_scan, plot_continuation

export FHNNetwork, complete_adjacency, ring_adjacency
export network_rhs!, network_rhs, synchronized_state, distributed_initial_state
export initial_state_from_orbit, simulate_network, network_state_matrix
export synchronization_error, lyapunov_dispersion, phase_order_squared, phase_order
export network_metrics, lyapunov_derivative, lyapunov_upper_bound, is_mean_field
export mean_field_linear_blocks, transverse_floquet
export network_layout, plot_network_snapshot, plot_network_timeseries
export plot_network_phase, plot_network_metrics, plot_coupling_comparison
export animate_network_dashboard
export save_network_result, load_network_result

# Etapa D — leyes suficientes, certificación y evaluación (sin optimizar J).
include("control/sufficient_controls.jl")
include("analysis/control_guarantees.jl")
include("analysis/control_costs.jl")
include("visualization/control_plots.jl")
include("utils/control_storage.jl")

export AttractiveModulatedCoupling
export FHNCurrentFeedback, feedback_currents, feedback_rhs!, feedback_rhs
export simulate_feedback, balanced_current_profile, balanced_current_law
export control_current_matrix, total_current_effort
export sufficient_decay_rate, sufficient_arrival_time, lyapunov_certificate
export feedback_lyapunov_derivative, feedback_lyapunov_upper_bound
export first_tolerance_hit, evaluate_sync_cost, evaluate_desync_cost
export plot_lyapunov_certificate, plot_control_comparison, plot_control_inputs
export plot_sync_cost_traces, plot_desync_cost_traces
export save_control_result, load_control_result


# Etapa E - Pontryagin, optimización directa acotada y diagnósticos.
include("control/optimal_pmp.jl")
include("control/optimal_direct.jl")
include("analysis/optimal_diagnostics.jl")
include("visualization/optimal_plots.jl")
include("utils/optimal_storage.jl")

export SyncOCPSpec, DesyncOCPSpec
export sync_error_gradient, phase_order_squared_gradient
export sync_ocp_field, desync_ocp_field
export sync_running_cost, desync_running_cost
export sync_hamiltonian, desync_hamiltonian
export sync_adjoint, desync_adjoint, sync_control_pmp, desync_control_pmp
export project_balanced_box, piecewise_value, piecewise_currents
export optimize_sync_direct, optimize_desync_direct
export simulate_optimal_sync, simulate_optimal_desync
export pmp_diagnostics_sync, pmp_diagnostics_desync, shoot_sync_pmp
export plot_optimal_sync, plot_optimal_desync, plot_pmp_sync, plot_pmp_desync
export save_optimal_result, load_optimal_result

# Etapa F - aproximaciones HJB en dimensiones pequeñas y comparaciones.
include("control/hjb_sync_reduced.jl")
include("control/hjb_desync_small.jl")
include("visualization/hjb_plots.jl")
include("visualization/hjb_comparison_video.jl")
include("utils/hjb_storage.jl")

export ReducedSyncHJB, solve_reduced_sync_hjb, reduced_sync_state, reduced_sync_rhs
export reduced_sync_value, reduced_sync_feedback, rollout_reduced_sync
export SmallDesyncHJB, solve_small_desync_hjb, small_desync_feedback, rollout_small_desync
export plot_hjb_sync, plot_hjb_desync, plot_hjb_comparison
export animate_stage_E_comparison, animate_stage_E_desync_comparison, save_hjb_result

end # module
