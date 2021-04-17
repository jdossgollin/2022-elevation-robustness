using DrWatson
@quickactivate "2021-elevation-robustness"

using Unitful
using HouseElevation
using Plots, UnitfulRecipes

function plot_elevation_cost()

    areas = [3_000, 2_000, 1_000]u"ft^2"
    elevations = (0:0.1:14)u"ft"

    p = plot(
        legend = :bottomright,
        xlabel = "Height Increase",
        ylabel = "Elevation Cost (Thousand \$USD)",
        ymin = 0,
        palette = :tab10,
        xticks = 0:2:14,
        size = (600, 400),
    )
    for area in areas
        house = HouseStructure(area, 185_000, 5u"ft") # last 2 don't matter here
        plot!(
            p,
            elevations,
            [elevation_cost(house, Δh) for Δh in elevations] ./ 1_000,
            label = "$area house",
            linewidth = 4,
        )
    end
    return p
end
p = plot_elevation_cost()
savefig(plotsdir("elevation_cost.pdf"))
