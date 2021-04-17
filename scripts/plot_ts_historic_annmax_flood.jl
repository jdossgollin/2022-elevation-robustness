#=
This script creates a time series of annual maximum floods
=#
using DrWatson
@quickactivate "2021-elevation-robustness"

using NorfolkFloods, Plots, Unitful, ColorSchemes

function plot_historic_annmax_flood()

    colors = ColorSchemes.tab10

    annual = get_norfolk_annual()
    surge_ft = (annual.msl + annual.surge) .|> x -> ustrip(u"ft", x)

    p = plot(
        annual.year,
        surge_ft,
        label = "",
        ylabel = "Annual Maximum Storm Surge (ft)",
        marker = :circle,
        size = (500, 300),
        title = "NOAA $(annual.gage_id): $(annual.gage_name)",
        titlefontsize = 10,
        labelfontsize = 9,
    )

    for (year, name, Δx, Δy, is_tc) in zip(
        [1933, 2003, 2015, 2009, 1962],
        ["Chesapeake–Potomac", "Isabel", "Irene", "Nor'Ida", "Ash Wednesday"],
        [12.5, -7.5, 2.5, 0, 10],
        [0.75, 0.5, 0.5, 1, 0.75],
        [true, true, true, false, false],
    )
        yobs = ustrip(surge_ft[findfirst(annual.year .== year)])
        x0 = year + Δx
        y0 = yobs + Δy
        color = is_tc ? colors[2] : colors[1]
        plot!(
            [x0, year],
            [y0, yobs],
            arrow = :closed,
            color = color,
            linewidth = 1,
            label = "",
        )
        scatter!(
            p,
            [x0],
            [y0],
            markercolor = :white,
            label = false,
            markerstrokecolor = :white,
            markersize = 10,
        )
        annotate!(p, x0, y0, text(name, :center, 7, color = color))
    end
    return p
end
p = plot_historic_annmax_flood()
savefig(p, plotsdir("ts_historic_anmax_flood.pdf"))
