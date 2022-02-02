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
    )
    plot!(p, ustrip.(u"ft", clearance), dmg_europa; linewidth=2, label="Europa (Not Used)")
    savefig(p, plots_dir("cost-expected-damage-emulator.pdf"))
    return p
end

function plot_cost_upfront(
    house_floor_area::T1, house_value_usd::T2, Δh_consider::Vector{<:Unitful.Length}
) where {T1<:Unitful.Area,T2<:Real}
    elevation_cost_fn = HouseElevation.get_elevation_cost_function()
    cost = elevation_cost_fn.(Δh_consider, house_floor_area)
    p = plot(
        ustrip.(u"ft", Δh_consider),
        cost ./ house_value_usd;
        linewidth=2,
        xlabel=L"Height Increase $\Delta h$ [ft]",
        ylabel="Up-Front Cost [% House Value]",
        leftmargin=5mm,
        label=false,
        xticks=0:2:14,
    )
    savefig(p, plots_dir("cost-up-front.pdf"))
    return p
end

function make_cost_plots(
    house_floor_area::T1, house_value_usd::T2, Δh_consider::Vector{<:Unitful.Length}
) where {T1<:Unitful.Area,T2<:Real}
    plot_depth_damage()
    plot_cost_expected_damage()
    plot_cost_upfront(house_floor_area, house_value_usd, Δh_consider)
    return nothing
end