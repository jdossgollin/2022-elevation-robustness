using LaTeXStrings
using Plots
using Plots: mm
using Unitful

"""Plot the depth damage relationship"""
function plot_depth_damage()
    dd_hazus = HouseElevation.get_depth_damage(:hazus)
    dd_europa = HouseElevation.get_depth_damage(:europa)

    depths = range(-5, 30; length=250)u"ft"

    p = plot(;
        legend=:topleft,
        xlabel="Flood Depth [ft relative to bottom of house]",
        ylabel="Damage [% House Value]",
        yformatter=pct_formatter,
    )
    plot!(
        p,
        ustrip.(u"ft", depths),
        dd_hazus.(depths);
        label="HAZUS (Used)",
        color=colors[1],
        linewidth=2,
    )
    plot!(
        p,
        ustrip.(u"ft", depths),
        dd_europa.(depths);
        label="Europa (Not Used)",
        color=colors[2],
        linewidth=2,
    )

    savefig(p, plots_dir("cost-depth-damage.pdf"))
    return p
end

function plot_cost_expected_damage()
    clearance = range(0, 15; length=250)u"ft"
    hazus = HouseElevation.get_expected_annual_damage_emulator(:hazus; overwrite=false)
    dmg_hazus = hazus.(clearance)
    europa = HouseElevation.get_expected_annual_damage_emulator(:europa; overwrite=false)
    dmg_europa = europa.(clearance)
    p = plot(
        ustrip.(u"ft", clearance),
        dmg_hazus;
        linewidth=2,
        xlabel="House Height [ft above MSL]",
        ylabel="Expected Annual Damage [% House Value]",
        leftmargin=5mm,
        label="HAZUS (Used)",
        color=colors[1],
        yformatter=pct_formatter,
    )
    plot!(
        p,
        ustrip.(u"ft", clearance),
        dmg_europa;
        linewidth=2,
        label="Europa (Not Used)",
        color=colors[2],
    )
    savefig(p, plots_dir("cost-expected-damage-emulator.pdf"))
    return p
end

function plot_cost_upfront(
    house_floor_area::T1, house_value_usd::T2, Δh_consider::Vector{<:Unitful.Length}
) where {T1<:Unitful.Area,T2<:Real}

    # get the function describing what it costs to elvate
    elevation_cost_fn = HouseElevation.get_elevation_cost_function()

    idx = findall(Δh_consider .> 0.0u"ft")
    Δh = Δh_consider[idx]
    cost = elevation_cost_fn.(Δh, house_floor_area)
    cost_prop = cost ./ house_value_usd

    p = plot(
        ustrip.(u"ft", Δh),
        cost_prop;
        linewidth=2,
        xlabel=L"Height Increase $\Delta h$ [ft]",
        ylabel="Up-Front Cost [% House Value]",
        leftmargin=5mm,
        label=false,
        xticks=0:2:14,
        color=colors[1],
    )
    # add in-line label
    annotate!(p, [8], [0.675], text("Cost of Elevating", :center, 10; color=colors[1]))

    scatter!(
        p,
        [0],
        [elevation_cost_fn(0.0u"ft", house_floor_area)];
        color=colors[2],
        label=false,
    )
    annotate!(
        p, [2.5], [0.05], text(L"No Cost if $\Delta h = 0$", :center, 10; color=colors[2])
    )

    savefig(p, plots_dir("cost-up-front.pdf"))
    return p
end
