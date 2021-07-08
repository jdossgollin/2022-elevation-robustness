### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 61c22ae1-ec2d-480c-b515-9cff34f3c4c3
begin
	using ColorSchemes
	using DataFrames
	using Distributions
	using DrWatson
	using Plots
	using StatsPlots
	using Turing: Chains
	using Unitful

	using Plots: mm
	using Printf: @sprintf
	
	using NorfolkBRICK
	using NorfolkFloods
end

# ╔═╡ 46688379-adb6-4ea7-8b2e-7d0864c652fc
md"""
load some required packages and set up the table of contents
"""

# ╔═╡ f0c4c942-ccb2-4776-a958-53bca97fc15f
md"""
# Observations

In this notebook we'll create a plot showing:

1. Historic floods
1. Return period estimates of storm surge
1. Projections of mean sea level (MSL) from the BRICK model (Wong et al, 2017) using both fast and slow ice sheet dynamics

> Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., & Keller, K. (2017). BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections. Geoscientific Model Development, 10(7), 2741–2760. https://doi.org/10.5194/gmd-10-2741-2017

"""

# ╔═╡ 66ac51c6-78ca-41bc-83ac-322867bec456
obs = NorfolkFloods.get_norfolk_annual(); # observed MSL and surge

# ╔═╡ 50174bdc-9c64-446a-a6b4-2866af053630
md"now define some constants"

# ╔═╡ b2864283-ce0b-45a5-9614-1fbbf9dcab26
colors = ColorSchemes.okabe_ito # colorblind friendly

# ╔═╡ 5e57a1a2-d4bc-4188-a4be-d47d1a39cddb
begin
	SYEAR = maximum(obs.year) + 1
	EYEAR = SYEAR + 100	
end;

# ╔═╡ 5e0fff62-acec-4b0f-860e-6fec91abb8f1
struct HistoricSurge
    name::String
    year::Int
    is_tc::Bool
end

# ╔═╡ 2236c1aa-ebf6-4138-ae8b-fa315d11938e
md"""
## Historic Floods

First, we will plot historic floods.
This combines both storm surge and mean sea level (MSL).
Notable events are demarcated on the plot!
"""

# ╔═╡ 8c9f8873-8183-4908-9b18-990f098f9475
function plot_historic_annmax_flood()

    annual = get_norfolk_annual()
    surge_ft = ustrip.(u"ft", annual.surge)

    p = plot(
        annual.year,
        surge_ft,
        label = "",
        ylabel = "Annual-Maximum Storm Surge (ft)",
        marker = :circle,
        markercolor = :gray,
        markerstrokewidth = 0,
        titlefontsize = 10,
        labelfontsize = 9,
    )

    storms = [
        HistoricSurge("Chesapeake–Potomac", 1933, true),
        HistoricSurge("Isabel", 2003, true),
        HistoricSurge("Irene", 2015, true),
        HistoricSurge("Nor'Ida", 2009, false),
        HistoricSurge("Ash Wednesday", 1962, false),
    ]

    for (storm, Δx, Δy) in zip(storms, [12.5, -7.5, 2.5, 0, 10], [0.75, 0.5, 1, 1, 0.75])
        yobs = ustrip(surge_ft[findfirst(annual.year .== storm.year)])
        x0 = storm.year + Δx
        y0 = yobs + Δy
        color = storm.is_tc ? colors[5] : colors[6]
        plot!(
            [x0, storm.year],
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
        annotate!(p, x0, y0, text(storm.name, :center, 7, color = color))
    end
    return p
end;

# ╔═╡ fdd9b856-de34-4280-8a28-804e8d583341
begin
	p1 = plot_historic_annmax_flood()
	title!(p1, "(a) MSL + Surge: Historic Floods")
	savefig(p1, plotsdir("historic_events.pdf"))
	plot(p1)
end

# ╔═╡ 537eac0e-3db8-4d0a-8771-3f2a1b035dbd
lsl_1933 = ustrip(u"ft", obs.msl[end] - obs.msl[findfirst(obs.year .== 1933)]);

# ╔═╡ 81d59382-3316-449c-99f7-0132f13ca6ab
md"""
The most extreme flood event on record is the Chesapeake-Potomac TC of 1933.
Since then there has been $(round(lsl_1933; digits=1)) ft of sea level rise.
If the same storm were to occur in $(obs.year[end]), then it would likely reach approximately $(round(lsl_1933 + ustrip(u"ft", obs.surge[findfirst(obs.year .== 1933)]); digits=1)) ft.
"""

# ╔═╡ e6d99cf8-479a-4001-b5d4-5c833c1554db
extract = function (traj::BRICKtrajectory)
    idx_f = findfirst(traj.years .== EYEAR)
    DataFrame(
        :lsl_m => ustrip(u"m", traj.lsl_m[idx_f]),
        :rcp => traj.rcp,
        :dynamics => traj.dynamics,
    )
end;

# ╔═╡ 7c23f843-8f05-4f86-87c5-ffda327d5386
md"""
## Return Period

Next we plot the return period.
The data is plotted using the empirical [Weibull Plotting Position](# https://glossary.ametsoc.org/wiki/Weibull_plotting_position) and the fit shown is the Bayesian storm surge model.
This model is discussed at great length in Notebook 02.
"""

# ╔═╡ 8fb3d1c4-7d9d-42e0-960f-967fd37e5795
function plotting_position(y::Vector{<:Real})
    N = length(y)
    rank = [sum(y .> yi) for yi in y] .+ 1
    aep = rank ./ (N + 1)
	rt = 1 ./ (aep)
	return rt
end;

# ╔═╡ 0b824539-24bf-4118-9614-8e518feb17ea
function plot_return_period()

    obs = get_norfolk_annual()
    y = ustrip.(u"ft", obs.surge)
    model = NorfolkFloods.GEVModel(y)
    posterior = get_fits(model, "stationary", 100_000, n_chains = 4)

    yhat = vcat(sample_predictive(posterior, length(obs.surge))[:]...)
    pp = plotting_position(y)

    rt_plot = 10 .^ (range(0, log10(125); length = 101)[2:end])
    quantile_plot = 1 .- 1 ./ rt_plot
    rl_plot = [quantile(yhat, pr) for pr in quantile_plot]
    rl_100 = quantile(yhat, 1 - 1 / 100)

    best_fits = hcat(
        [
            quantile(GeneralizedExtremeValue(μ, σ, ξ), quantile_plot) for
            (μ, σ, ξ) in zip(posterior[:μ], posterior[:σ], posterior[:ξ])
        ]...,
    )

    xticks = [2, 5, 10, 25, 50, 100]
    ticklabels = [@sprintf("%d", x) for x in xticks]
    p = plot(
        rt_plot,
        rl_plot,
        formatter = identity,
        xticks = (xticks, ticklabels),
        xaxis = :log,
        legend = :topleft,
        label = "Best Fit",
        xlabel = "Return Period (years)",
        ylabel = "Storm Surge (ft)",
        titlefontsize = 10,
        labelfontsize = 9,
    )
    for q in ((0.1, 0.9), (0.025, 0.975))
        lower = [quantile(best_fits[i, :], minimum(q)) for i = 1:size(best_fits)[1]]
        upper = [quantile(best_fits[i, :], maximum(q)) for i = 1:size(best_fits)[1]]
        ci = Int(round(100 * (maximum(q) - minimum(q))))
        plot!(
            p,
            rt_plot,
            lower,
            fillrange = upper,
            fillalpha = 0.3,
            label = "$(ci)% CI",
            color = "gray",
            linewidth = 0,
        )
    end
    scatter!(
        p,
        pp,
        y,
        color = :black,
        label = "Observed",
        markerstrokewidth = 0,
        markersize = 3,
    )
    Plots.vline!(p, [100], linestyle = :dash, color = :gray, label = false)
    return p
end

# ╔═╡ b9e764d0-a2c7-492b-add3-e8d480c1643b
begin
	p2 = plot_return_period()
    title!(p2, "(b) Surge Only: Stationary GEV Return Periods")
	savefig(p2, plotsdir("surge_return_period.pdf"))
	plot(p2)
end

# ╔═╡ 8f23aceb-d69a-4b9b-9420-bd347fc797b9
md"""
According to this Bayesian GEV analysis, the storm surge associated with the Chesapeake-Potomac TC had a return period just over 100 years.
"""

# ╔═╡ a097e3c0-06cb-43a3-aadc-00b87094c1a4
md"""
## Sea Level Rise

Finally we plot projections of mean sea level
"""

# ╔═╡ 20f74cc5-99e7-49db-8dbc-4f2744018d1e
function plot_msl(dynamics::String; q = [0.95, 0.05])

    obs = get_norfolk_annual()
    lsl_trajectories = get_norfolk_brick(syear = minimum(obs.year); eyear = EYEAR)

    p = plot(
        xlabel = "Year",
        ylabel = "MSL at Norfolk, VA (ft)",
        legend = :topleft,
        titlefontsize = 10,
        labelfontsize = 9,
    )
    scatter!(
        p,
        obs.year,
        ustrip.(u"ft", obs.msl),
        label = "Observed",
        color = :black,
        markerstrokewidth = 0,
        markersize = 3,
    )
    for (rcp, color) in zip([2.6, 4.5, 6.0, 8.5], colors)
        trajs = [
            traj for
            traj in lsl_trajectories if (traj.rcp == rcp) & (traj.dynamics == dynamics)
        ]
        sims = ustrip.(u"ft", hcat([traj.lsl_m for traj in trajs]...))
        N = size(sims)[1]
        upper = [quantile(sims[i, :], maximum(q)) for i = 1:N]
        lower = [quantile(sims[i, :], minimum(q)) for i = 1:N]
        plot!(
            p,
            trajs[1].years,
            lower,
            fillrange = upper,
            fillalpha = 0.1,
            label = "RCP $rcp",
            color = color,
            linestyle = :solid,
            linewidth = 2,
        )
        plot!(
            p,
            trajs[1].years,
            upper,
            label = false,
            color = color,
            linewidth = 2,
            linestyle = :dash,
        )
    end
    return p
end;

# ╔═╡ 0b88699a-28de-464b-98cb-3dc5b9ee865f
begin
	p3 = plot_msl("fast")
    title!(p3, "(c) MSL Only: Fast Ice Sheet Dynamics")
	savefig(p3, plotsdir("msl_brick_fast.pdf"))
	plot(p3)
end

# ╔═╡ 97935eb9-a9e5-4df0-b968-2e97c2393efc
md"""
There are a lot of lines on this figure.

1. Solid lines are lower bounds of the 90% confidence interval
1. Dashed lines are upper bounds of the 90% confidence interval
1. Colors correspond to 4 RCP scenarios
"""

# ╔═╡ 8fbea3b1-ddf6-4ba6-ad2e-82754d88cd66
begin
	p4 = plot_msl("slow")
    title!(p4, "(d) MSL Only: Slow Ice Sheet Dynamics")
    plot!(p4, legend = false) # only need one
	savefig(p4, plotsdir("msl_brick_slow.pdf"))
	plot(p4)
end

# ╔═╡ 999e51f5-ccb4-42d7-a00c-2a008eb75f9e
md"""
We can see that under both fast and slow ice sheet dynamics, mean sea levels are relatively stable through about 2050; different scenarios diverge rapidly in the mid to late 21st century.
"""

# ╔═╡ 821b8cd7-540e-4f96-814d-281a7308e308
md"""
Now let's put those figures together to create a final figure for sharing.
"""

# ╔═╡ 2f869544-f744-4082-b883-6be1ea5e3267
begin
	ylims!(p1, (2.5, 8.5))
	ylims!(p2, (2.5, 8.5))
	ylabel!(p2, "")
	ylabel!(p4, "")
	p = plot(
		p1,
		p2,
		p3,
		p4;
		size = (1000, 750),
		bottom_margin = 5mm,
		left_margin = 5mm,
		top_margin = 2.5mm,
		link = :y,
	)
end

# ╔═╡ ffdfc2e8-6b1a-452b-9f7a-75479663840d
savefig(p, plotsdir("historic_quad.pdf"))

# ╔═╡ Cell order:
# ╟─46688379-adb6-4ea7-8b2e-7d0864c652fc
# ╠═61c22ae1-ec2d-480c-b515-9cff34f3c4c3
# ╟─f0c4c942-ccb2-4776-a958-53bca97fc15f
# ╠═66ac51c6-78ca-41bc-83ac-322867bec456
# ╟─50174bdc-9c64-446a-a6b4-2866af053630
# ╠═b2864283-ce0b-45a5-9614-1fbbf9dcab26
# ╠═5e57a1a2-d4bc-4188-a4be-d47d1a39cddb
# ╠═5e0fff62-acec-4b0f-860e-6fec91abb8f1
# ╟─2236c1aa-ebf6-4138-ae8b-fa315d11938e
# ╠═8c9f8873-8183-4908-9b18-990f098f9475
# ╠═fdd9b856-de34-4280-8a28-804e8d583341
# ╠═537eac0e-3db8-4d0a-8771-3f2a1b035dbd
# ╟─81d59382-3316-449c-99f7-0132f13ca6ab
# ╠═e6d99cf8-479a-4001-b5d4-5c833c1554db
# ╟─7c23f843-8f05-4f86-87c5-ffda327d5386
# ╠═8fb3d1c4-7d9d-42e0-960f-967fd37e5795
# ╠═0b824539-24bf-4118-9614-8e518feb17ea
# ╠═b9e764d0-a2c7-492b-add3-e8d480c1643b
# ╟─8f23aceb-d69a-4b9b-9420-bd347fc797b9
# ╟─a097e3c0-06cb-43a3-aadc-00b87094c1a4
# ╠═20f74cc5-99e7-49db-8dbc-4f2744018d1e
# ╠═0b88699a-28de-464b-98cb-3dc5b9ee865f
# ╟─97935eb9-a9e5-4df0-b968-2e97c2393efc
# ╠═8fbea3b1-ddf6-4ba6-ad2e-82754d88cd66
# ╟─999e51f5-ccb4-42d7-a00c-2a008eb75f9e
# ╟─821b8cd7-540e-4f96-814d-281a7308e308
# ╠═2f869544-f744-4082-b883-6be1ea5e3267
# ╠═ffdfc2e8-6b1a-452b-9f7a-75479663840d
