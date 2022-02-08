using HouseElevation
using Plots
using Plots: mm
using StatsPlots

"""Plot some time-varying PDFs"""
function plot_lsl_pdfs(
    all_trajs::Vector{<:HouseElevation.BRICKSimulation}; ub=0.95, lb=0.05
)
    width = ub - lb
    p = plot(;
        legend=:left,
        xlabel="Time [year]",
        ylabel="Mean Sea Level at Sewells Point, VA [ft]",
        legendtitle="$(Int(ceil(width * 100)))% Interval:",
    )

    # which PDFs to plot?
    models = ((2.6, "slow"), (4.5, "fast"), (8.5, "fast"))

    for (model, color) in zip(models, colors[3:(length(models) + 2)])
        (rcp, dyn) = model
        trajs = [traj for traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == dyn)]

        years = first(trajs).years

        sims = ustrip.(u"ft", hcat([traj.lsl for traj in trajs]...))
        N = size(sims)[1]

        upper = [quantile(sims[i, :], ub) for i in 1:N]
        lower = [quantile(sims[i, :], lb) for i in 1:N]
        plot!(
            p,
            years,
            upper;
            fillrange=lower,
            fillalpha=0.5,
            label="RCP $rcp / $dyn BRICK",
            color=color,
            linewidth=0,
            leftmargin=5mm,
        )
    end
    vline!(p, [2100]; linestyle=:dash, label=false, color=:gray, linewidth=2)
    return p
end

"""Plot some booxplots"""
function plot_lsl_boxplots_2100(all_trajs::Vector{<:HouseElevation.BRICKSimulation})
    all_rcp = unique([traj.rcp for traj in all_trajs])

    p = plot(; xticks=(1:length(all_rcp), "RCP " .* string.(all_rcp)), legend=:bottomright)
    for (rcp_idx, rcp) in enumerate(all_rcp)
        msl_slow = [
            get_year_data(traj, 2100) for
            traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == "slow")
        ]
        msl_fast = [
            get_year_data(traj, 2100) for
            traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == "fast")
        ]
        if rcp_idx == 1
            slow_label = "Slow BRICK Dynamics"
            fast_label = "Fast BRICK Dynamics"
        else
            slow_label = false
            fast_label = false
        end
        violin!(
            p,
            [rcp_idx],
            ustrip.(u"ft", msl_slow);
            color=colors[1],
            side=:left,
            label=slow_label,
        )
        violin!(
            p,
            [rcp_idx],
            ustrip.(u"ft", msl_fast);
            color=colors[2],
            side=:right,
            label=fast_label,
        )
    end
    return p
end

"""Make all the LSL plots"""
function plot_brick(all_trajs::Vector{<:HouseElevation.BRICKSimulation})
    p1 = plot_lsl_pdfs(all_trajs)
    p2 = plot_lsl_boxplots_2100(all_trajs)
    add_panel_letters!([p1, p2]; fontsize=14)
    p = plot(
        p1,
        p2;
        link=:y,
        ylims=(0, 6.5),
        layout=grid(1, 2; widths=[0.7, 0.3]),
        size=[1000, 400] .* 1.125,
        xtickfontsize=9,
        ylabelfontsize=9,
        bottom_margin=5mm,
        leftmargin=5mm,
    )
    savefig(p, plots_dir("lsl-evolution.pdf"))
    return p
end
