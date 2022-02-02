using Revise
using HouseElevation

using ColorSchemes
using Unitful

# HouseElevation.clear_cache() # to clean all the processed data and plots
colors = ColorSchemes.okabe_ito # colorblind friendly and consistent scheme

include("01-surge-modeling.jl")
include("02-mean-sea-level.jl")
include("03-cost.jl")

function main()
    start_year = 2022
    end_year = 2071

    # make plots of the annual surge data
    stn = HouseElevation.TidesAndCurrentsRecord()
    annual = HouseElevation.get_annual_data(stn)
    fits = get_norfolk_posterior()
    make_surge_plots(annual, fits)

    # make some plots of the sea level (BRICK) data
    all_trajs = HouseElevation.get_norfolk_brick(; syear=2022, eyear=2100)
    make_lsl_plots(all_trajs)

    # make some plots of the cost functions
    house_floor_area = 1500u"ft^2"
    house_value_usd = 200_000
    elevation_init = 7u"ft"
    discount_rate = 0.02
    Δh_consider = collect(0:0.25:14)u"ft" # this is the decision space!
    make_cost_plots(house_floor_area, house_value_usd, Δh_consider)

    # almost done -- just need to clean this bit up
    s = HouseElevation.get_norfolk_brick(; syear=start_year, eyear=end_year)
    f = HouseElevation.get_system_model(s)
    u = HouseElevation.exhaustive_exploration(f, s, x)
    HouseElevation.total_cost_usd.(u)

    return nothing
end

#main()
#supplemental()
