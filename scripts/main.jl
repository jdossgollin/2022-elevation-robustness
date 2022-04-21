using Revise
using HouseElevation

using ColorSchemes
using Unitful

# HouseElevation.clear_cache() # to clean all the processed data and plots
include("plotutils.jl")
colors = get_colormap()

# these files just provide functions that will make plots -- nothing happens except defining functions
include("01-surge-modeling.jl")
include("02-mean-sea-level.jl")
include("03-cost.jl")
include("04-scenario-maps.jl")
include("05-tradeoffs-by-rcp.jl")
include("06-combining.jl")

function main()

    # define some constants
    house_floor_area = 1500u"ft^2"
    house_value_usd = 200_000.0 # HOUSE NOT LAND VALUE
    discount_rate = 0.025 # mortgages going for 3-4.5% at the moment
    syear = 2022
    eyear = 2092
    # for a **VERY** rough idea of prices see
    # https://www.zillow.com/homedetails/9638-Selby-Pl-Norfolk-VA-23503/79223088_zpid/

    # get the storm surge posterior samples
    fits = get_surge_posterior()

    # make plots of the annual surge data
    let
        stn = HouseElevation.TidesAndCurrentsRecord()
        annual = HouseElevation.get_annual_data(stn)
        plot_surge_gev_priors()
        plot_surge_prior_chains()
        plot_surge_synthetic_experiment(annual)
        plot_surge_posterior_chains(fits)
        plot_surge_posterior_teststats(annual, fits)
        plot_surge_posterior_return(annual, fits)
        plot_surge_obs_return(annual, fits)
        plot_surge_prior_return()
    end

    # make some plots of the sea level (BRICK) data
    # uses notation from paper: s = set of all scenarios
    let
        s = HouseElevation.get_lsl(; syear=syear, eyear=2125)
        plot_lsl_evolution(s)
    end

    # make some plots of the cost functions
    let
        Δh_consider = collect(0:0.25:14)u"ft" # this is the decision space!
        plot_depth_damage()
        plot_cost_expected_damage()
        plot_cost_upfront(house_floor_area, house_value_usd, Δh_consider)
    end

    # create scenario maps
    let
        elevation_init_bfe = [-2, -0.5, 1]u"ft" # height relative to gauge
        x_plot = [0, 4, 8]u"ft" # elevations to plot
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
    end

    let
        s = HouseElevation.get_lsl(; syear=syear, eyear=eyear)
        bfe = calc_bfe(fits, s, syear)
        x = collect(0:0.25:12)u"ft"
        elevation_init = bfe - 1u"ft"
        u = get_outcomes(
            x;
            syear=syear,
            eyear=eyear,
            fits=fits,
            house_floor_area=house_floor_area,
            elevation_init=elevation_init,
            discount_rate=discount_rate,
            house_value_usd=house_value_usd,
            overwrite=false,
        )
        plot_scenario_map_height_slr(; x=x, s=s, u=u, house_value_usd=house_value_usd)

        # plot tradeoffs by RCP scenario
        plot_rcp_tradeoffs(
            u, s, x; house_value_usd=house_value_usd, house_floor_area=house_floor_area
        )

        # write the priors to a table
        HouseElevation.write_priors(plots_dir("lsl_priors.tex"))

        # plot the implicit weight
        plot_priors_weights(s)

        # plot tradeoffs by prior
        plot_prior_tradeoffs(
            u, s, x; house_value_usd=house_value_usd, house_floor_area=house_floor_area
        )
    end

    plot_grid_scheme()

    return nothing
end

main()
