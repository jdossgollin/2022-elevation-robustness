using Revise
using HouseElevation

using ColorSchemes

# HouseElevation.clear_cache() # to clean all the processed data and plots
colors = ColorSchemes.okabe_ito # colorblind friendly and consistent scheme

include("01-surge-modeling.jl")

function supplemental()

    # STORM SURGE MODELING
    plot_surge_gev_priors()
    plot_annmax_floods()
    plot_surge_prior_chains()
    plot_surge_synthetic_experiment()
    plot_surge_posterior_chains()
    plot_surge_posterior_teststats()
    plot_surge_posterior_return()
    return nothing
end

function main()
    return plot_fig1()
end

main()
# supplemental()