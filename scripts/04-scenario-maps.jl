using LaTeXStrings
using Plots
using Plots: mm
using Unitful

function plot_scenario_map_slr_cost(;
    u_plot::Matrix{<:HouseElevation.Outcome},
    x_plot::Vector{<:Unitful.Length},
    s_plot::Vector{<:HouseElevation.BRICKSimulation},
    house_value_usd,
)
    # get the amount of SLR over the time domain
    s1 = first(s_plot)
    syear = minimum(s1.years)
    eyear = maximum(s1.years)
    msl_rise = get_year_data(sows, eyear) .- get_year_data(sows, syear)

    plots = []
    for (i, Δh) in enumerate(x_plot)
        total_cost_prop =
            [outcome.led_usd + outcome.upfront_cost_usd for outcome in outcomes[:, i]] ./ house_value_usd .*
            100
        title = "($('a'+(i-1))): Δh = $Δh"
        label = i == 1 ? "1 dot = 1 SOW" : false
        p = scatter(
            ustrip.(u"ft", msl_rise),
            total_cost_prop;
            markersize=0.5,
            markercolor=:black,
            xlabel="Local Sea Level Rise from $syear to $eyear [ft]",
            label=label,
            title_align=:left,
            title=title,
            ylabel="Total Cost [% House Value]",
            legend=:bottomright,
        )
        push!(plots, p)
    end
    p = plot(plots...; size=(1000, 750), leftmargin=5mm, link=:y, ylims=(50, 350), dpi=300)
    savefig(p, plots_dir("scenario-map-slr-cost.png"))
    return p
end

function make_scenario_maps(;
    u_plot::Matrix{<:HouseElevation.Outcome},
    x_plot::Vector{<:Unitful.Length},
    s_plot::Vector{<:HouseElevation.BRICKSimulation},
    house_value_usd,
)
    plot_scenario_map_slr_cost(;
        u_plot=u_plot, x_plot=x_plot, s_plot=s_plot, house_value_usd=house_value_usd
    )
    return nothing
end
