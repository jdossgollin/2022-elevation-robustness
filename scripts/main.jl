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
    elevation_init = 7u"ft"
    discount_rate = 0.015
    start_year = 2022
    end_year = start_year + 70

    # make the scenario maps
    u_plot, s_plot, x_plot = get_outcomes(
        [0.0, 3, 6, 9]u"ft";
        syear=start_year,
        eyear=end_year,
        fits=fits,
        house_floor_area=house_floor_area,
        elevation_init=elevation_init,
        discount_rate=discount_rate,
        house_value_usd=house_value_usd,
        overwrite=false,
    )
    make_scenario_maps(;
        u_plot=u_plot, s_plot=s_plot, x_plot=x_plot, house_value_usd=house_value_usd
    )

    return nothing
end

main()