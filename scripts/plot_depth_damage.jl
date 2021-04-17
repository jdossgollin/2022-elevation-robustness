using DrWatson
@quickactivate "2021-elevation-robustness"

using Unitful
using HouseElevation
using Plots, UnitfulRecipes

function plot_depth_damage()

    depths = (-2:0.1:16)u"ft"
    house = HouseStructure(1000u"ft^2", 185_000, 0u"ft")

    p = plot(
        legend = :bottomright,
        xlabel = "Flood Depth Relative to House",
        ylabel = "Damages (Proportion of House Value)",
        ymin = 0,
        palette = :tab10,
        size = (600, 400),
    )
    for model in [:hazus, :europa]
        plot!(
            p,
            depths,
            [depth_damage_frac(house, depth, model) for depth in depths],
            label = string(model),
            linewidth = 4,
        )
    end
    return p
end

p = plot_depth_damage()
savefig(plotsdir("depth_damage.pdf"))