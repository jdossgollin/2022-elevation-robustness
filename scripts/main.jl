using Revise
using HouseElevation

using ColorSchemes

# HouseElevation.clear_cache() # to clean all the processed data and plots
colors = ColorSchemes.okabe_ito # colorblind friendly and consistent scheme

include("01-surge-modeling.jl")
include("02-mean-sea-level.jl")
include("03-cost.jl")

function supplemental()

    # STORM SURGE MODELING
    plot_surge_gev_priors()
    plot_annmax_floods()
    plot_surge_prior_chains()
    plot_surge_synthetic_experiment()
    plot_surge_posterior_chains()
    plot_surge_posterior_teststats()
    plot_surge_posterior_return()
    # COSTS
    plot_depth_damage()
    plot_cost_expected_damage()
    return nothing
end

function main()
    plot_surge_obs_return()
    plot_lsl_pdfs()
    plot_lsl_boxplots_2100()
    return nothing
end

main()
#supplemental()