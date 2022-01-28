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
    model = HouseElevation.get_expected_damage_emulator(:hazus)
    clearance = range(0, 15; length=250)u"ft"
    dmg = model.(clearance)
    p = plot(
        ustrip.(u"ft", clearance),
        dmg;
        linewidth=2,
        xlabel="House Height [ft above MSL]",
        ylabel="Expected Annual Damage [% House Value]",
        leftmargin=5mm,
        label=false,
    )
    savefig(p, plots_dir("cost-expected-damage-emulator.pdf"))
    return p
end
