module FitzHughNagumoControl

using LinearAlgebra
using OrdinaryDiffEq
using CairoMakie
using JLD2
using Dates

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

end # module
