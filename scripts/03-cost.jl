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
    hazus = HouseElevation.get_expected_damage_emulator(:hazus; overwrite=false)
    dmg_hazus = hazus.(clearance)
    europa = HouseElevation.get_expected_damage_emulator(:europa; overwrite=false)
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

function plot_cost_upfront()
    A = 1_000u"ft^2"
    value = 185_000 # house value in dollars -- see Z21
    Δhs = range(0, 14; length=100)u"ft"
    cost = HouseElevation.elevation_cost.(A, Δhs)
    p = plot(
        ustrip.(u"ft", Δhs),
        cost ./ value;
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
