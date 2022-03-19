using HouseElevation
using ColorSchemes
using DataFrames
using Plots
using Plots: mm
using StatsPlots

"""Plot some time-varying PDFs"""
function plot_lsl_pdfs(all_trajs::Vector{<:HouseElevation.LSLSim}; ub=0.95, lb=0.05)
    width = ub - lb
    p = plot(;
        xlabel="Time [year]",
        ylabel="Mean Sea Level at Sewells Point, VA [ft]",
        legendtitle="$(Int(ceil(width * 100)))% Interval:",
        legend=(0.2, 0.8),
    )

    # which PDFs to plot?
    models = (
        (rcp=2.6, model="BRICK Slow"), (rcp=4.5, model="K14"), (rcp=8.5, model="DP16")
    )

    for (model, color) in zip(models, ColorSchemes.seaborn_bright6)
        (rcp, modelname) = model
        trajs = [
            traj for traj in all_trajs if (traj.rcp == rcp) & (traj.model == modelname)
        ]

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
            label="RCP $rcp, $modelname",
            color=color,
            linewidth=0,
            leftmargin=5mm,
        )
    end
    vline!(p, [2100]; linestyle=:dash, label=false, color=:gray, linewidth=2)
    return p
end

"""Plot some booxplots"""
function plot_lsl_boxplots_2100(all_trajs::Vector{<:HouseElevation.LSLSim})
    msl_2100 = vcat(
        [
            DataFrame(
                "RCP" => traj.rcp,
                "Model" => traj.model,
                "msl_2100_ft" => ustrip(u"ft", get_year_data(traj, 2100)),
            ) for traj in all_trajs
        ]...,
    )
    p = @df msl_2100 groupedviolin(
        :Model,
        :msl_2100_ft,
        group=:RCP,
        palette=colors,
        legend=(0.15, 0.8),
        legendtitle="  RCP  ",
        ylabel="Mean Sea Level in 2100",
        xrotation=30,
    )
    return p
end

"""Make all the LSL plots"""
function plot_lsl_evolution(all_trajs::Vector{<:HouseElevation.LSLSim})
    p1 = plot_lsl_pdfs(all_trajs)
    p2 = plot_lsl_boxplots_2100(all_trajs)
    add_panel_letters!([p1, p2]; fontsize=14)
    p = plot(
        p1,
        p2;
        link=:y,
        ylims=(-0.5, 12),
        layout=grid(1, 2; widths=[0.6, 0.4]),
        size=[1000, 400] .* 1.125,
        xtickfontsize=9,
        ylabelfontsize=9,
        bottom_margin=7.5mm,
    )
    savefig(p, plots_dir("lsl-evolution.pdf"))
    return p
end
