using DrWatson
@quickactivate "2021-elevation-robustness"

using NorfolkFloods, Plots, Unitful

annual = get_norfolk_annual()
surge_ft = (annual.msl + annual.surge) .|> x -> uconvert(u"ft", x) .|> ustrip
plot(
    annual.year,
    surge_ft,
    label = false,
    xlabel = "Year",
    ylabel = "Annual Maximum Storm Surge (ft)",
    size = (600, 450),
    title = "NOAA $(annual.gage_id): $(annual.gage_name)",
    marker = :circle,
)
annotate!(1944, 7, text("Chesapeake–Potomac", :darkorange2, :left, 7))
plot!(
    [1943.5, 1933],
    [7, ustrip(surge_ft[findfirst(annual.year .== 1933)])],
    arrow = true,
    color = :darkorange2,
    linewidth = 1,
    label = "",
)

annotate!(1998, 6.75, text("Isabel", :darkorange2, :right, 7))
plot!(
    [1998.25, 2003],
    [6.75, ustrip(surge_ft[findfirst(annual.year .== 2003)])],
    arrow = true,
    color = :darkorange2,
    linewidth = 1,
    label = "",
)

annotate!(2015, 6.75, text("Irene", :darkorange2, :center, 7))
plot!(
    [2015, 2011],
    [6.65, ustrip(surge_ft[findfirst(annual.year .== 2011)])],
    arrow = true,
    color = :darkorange2,
    linewidth = 1,
    label = "",
)

annotate!(2009, 7.1, text("Nor'Ida", :purple, :center, 7))
plot!(
    [2009, 2009],
    [7, ustrip(surge_ft[findfirst(annual.year .== 2009)])],
    arrow = true,
    color = :purple,
    linewidth = 1,
    label = "",
)

annotate!(1968, 6.35, text("Ash Wednesday\nNor'Easter", :purple, :left, 7))
plot!(
    [1967.8, 1962],
    [6.35, ustrip(surge_ft[findfirst(annual.year .== 1962)])],
    arrow = true,
    color = :purple,
    linewidth = 1,
    label = "",
)

savefig(plotsdir("ann_max_floods.pdf"))
