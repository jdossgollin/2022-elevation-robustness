using Revise
using HouseElevation

using ColorSchemes
using Unitful

# HouseElevation.clear_cache() # to clean all the processed data and plots
colors = ColorSchemes.okabe_ito # colorblind friendly and consistent scheme

include("plotutils.jl")
include("01-surge-modeling.jl")
include("02-mean-sea-level.jl")
include("03-cost.jl")
include("04-scenario-maps.jl")

function main()

    # make plots of the annual surge data
    stn = HouseElevation.TidesAndCurrentsRecord()
    annual = HouseElevation.get_annual_data(stn)
    fits = get_norfolk_posterior()
    make_surge_plots(annual, fits)

    # make some plots of the sea level (BRICK) data
    all_trajs = HouseElevation.get_norfolk_brick(; syear=2022, eyear=2122)
    make_lsl_plots(all_trajs)

    # make some plots of the cost functions
    house_floor_area = 1500u"ft^2"
    house_value_usd = 200_000.0
    Δh_consider = collect(0:0.25:14)u"ft" # this is the decision space!
    make_cost_plots(house_floor_area, house_value_usd, Δh_consider)

    # some more constants
    discount_rate = 0.015
    syear = 2022
    eyear = syear + 70

    # make the scenario maps
    elevation_init_bfe = [-5, -2, 1]u"ft" # height relative to gauge
    x_plot = [0, 3, 6, 9]u"ft"

    plot_scenario_map_slr_cost(;
        x_plot=x_plot,
        elevation_init_bfe=elevation_init_bfe,
        fits=fits,
        house_value_usd=house_value_usd,
        syear=syear,
        eyear=eyear,
        house_floor_area=house_floor_area,
        discount_rate=discount_rate,
        overwrite=false,
    )

    return nothing
end

main()
