using Revise
using HouseElevation

using ColorSchemes
using Unitful

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
    plot_cost_upfront()
    return nothing
end

function main()
    plot_surge_obs_return()
    plot_lsl_pdfs()
    plot_lsl_boxplots_2100()
    return nothing
end

#main()
#supplemental()

start_year = 2022
end_year = 2071
s = HouseElevation.get_norfolk_brick(; syear=start_year, eyear=end_year)
f = HouseElevation.get_system_model(s)
x = collect(0:2:14)u"ft"
si = first(sows)
xj = rand(x)
u_ij = f(si, xj)

u = HouseElevation.exhaustive_exploration(f, s, x)
u[1, 1]