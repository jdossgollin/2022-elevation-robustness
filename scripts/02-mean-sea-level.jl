using HouseElevation
using Plots
using Plots: mm
using StatsPlots

"""Plot some time-varying PDFs"""
function plot_lsl_pdfs(
    all_trajs::Vector{<:HouseElevation.BRICKSimulation}; ub=0.95, lb=0.05
)
    p = plot(;
        legend=:topleft,
        xlabel="Time [year]",
        ylabel="Mean Sea Level at Sewells Point, VA [ft]",
    )

    # which PDFs to plot?
    models = ((2.6, "slow"), (6.0, "fast"), (8.5, "fast"))

    for (model, color) in zip(models, colors)
        rcp, dyn = model
        trajs = [traj for traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == dyn)]

        years = first(trajs).years

        sims = ustrip.(u"ft", hcat([traj.lsl for traj in trajs]...))
        N = size(sims)[1]

        upper = [quantile(sims[i, :], ub) for i in 1:N]
        lower = [quantile(sims[i, :], lb) for i in 1:N]
        width = ub - lb
        plot!(
            p,
            years,
            upper;
            fillrange=lower,
            fillalpha=0.5,
            label="RCP $rcp / $dyn BRICK, $(Int(ceil(width * 100)))% CI",
            color=color,
            linewidth=0,
            leftmargin=5mm,
        )
    end
    savefig(p, plots_dir("lsl-pdfs-timevarying.pdf"))
    return p
end

"""Plot some booxplots"""
function plot_lsl_boxplots_2100(all_trajs::Vector{<:HouseElevation.BRICKSimulation})
    df = vcat(
        [
            DataFrame(
                "scenario" => "RCP $(traj.rcp) / $(traj.dynamics) BRICK",
                "msl_2100" => ustrip(u"ft", get_year_data(traj, 2100)),
            ) for traj in all_trajs
        ]...,
    )
    p = plot(;
        ylabel="MSL at Sewells Point, VA in 2100 [ft]",
        bottom_margin=10mm,
        leftmargin=10mm,
        left_margin=5mm,
        size=(700, 400),
        xtickfontsize=9,
        ylabelfontsize=9,
    )
    @df df violin!(p, :scenario, :msl_2100, xrotation=35, legend=false)
    savefig(p, plots_dir("lsl-boxplot-2100.pdf"))
    return p
end

"""Make all the LSL plots"""
function make_lsl_plots(all_trajs::Vector{<:HouseElevation.BRICKSimulation})
    plot_lsl_pdfs(all_trajs)
    plot_lsl_boxplots_2100(all_trajs)
    return nothing
end
