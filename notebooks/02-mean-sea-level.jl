### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 6eaeeaaf-5fe8-44b2-8957-b47918f26b17
begin
	# we start by loading in the packages we use
	using Pkg
	Pkg.activate("..")
	using ColorSchemes
	using CSV
	using DataFrames
	using DrWatson
	using NetCDF
	using PlutoUI
	using Plots
	using Plots: mm
	using Statistics
	using StatsBase
	using StatsPlots
	using Unitful
	using UnitfulRecipes
	TableOfContents()
end

# ╔═╡ 49cd8998-378d-11ec-3ceb-811bf9c58697
md"""
# Notebook 02: Mean Sea Level

> Doss-Gollin and Keller, 2021, In Prep. **Do not circulate!**

In this notebook we'll explore simulations of mean sea level.

*If you have not already run Notebook 01, stop and do that!*
"""

# ╔═╡ ab567401-943f-4877-8572-cfa8ffb531d4
md"""
## Read Data
"""

# ╔═╡ c259b7d3-0873-4502-9d4b-a5d18c5adcd8
md"""
The main projections of sea level rise that we use come from the BRICK model with two parameterizations of ice sheet dynamics ('fast' and 'slow') and four RCP scenarios.
For a detailed comparison of these simulations to other model structures and assumptions, see

> Ruckert, K. L., Srikrishnan, V., & Keller, K. (2019). Characterizing the deep uncertainties surrounding coastal flood hazard projections: A case study for Norfolk, VA. Scientific Reports, 9(1), 1–12. https://doi.org/10.1038/s41598-019-47587-6
"""

# ╔═╡ 5048e15a-d382-4a3c-9aa6-d0d34cb2b297
md"""
## Plots

Let's get straight into the analysis
"""

# ╔═╡ a4a6d097-0d5c-48c7-b8bc-34baa3c910cc
colors = ColorSchemes.okabe_ito; # colorblind friendly

# ╔═╡ f33c8247-c110-412d-a7d4-4b8aa8f89f34
md"""
First, let's illustrate the 'multiple PDF problem' by plotting all the PDFs for the year 2100
"""

# ╔═╡ 37d98257-5de6-456a-a542-36d0dc1c7a96
md"""
In this (admittedly noisy!) plot we can see different CDFs of sea level in 2100.
We see that the ice sheet dynamics ('fast' vs 'slow') don't matter much (prior to 2100) for low emissions scenarios, but they make a big difference for the higher emissions pathway and at the tails of the distribution.
"""

# ╔═╡ f251e057-4409-4c8b-a2d7-4603c8af6750
md"""
Next, let's plot 95% Confidence Intervals for a couple of representative models
"""

# ╔═╡ 5f912aeb-c9e4-485a-897f-a71e84ec6038
md"""
Alternatively, we can show individual scenarios with lines
"""

# ╔═╡ f815e33f-97fa-4a1b-b8fc-77241d81d172
md"""
Last, we can create boxplots
"""

# ╔═╡ 48b37aa3-b756-4628-be82-5289ed202aa2
md"""
## Save Simulations

We'll save our trajectories for easy use later
"""

# ╔═╡ 15ad5394-ee85-4aff-851d-e9e9a58ff1a0
md"""
## Helper Functions

The functions below help us work with the code we need.
Pluto notebooks do some fancy tricks to let us store this code at the bottom of the file, even though we need to run them before we run some of the other blocks.
"""

# ╔═╡ 3a5a60eb-3493-4359-924e-ad19ddd28b63
md"""
A single BRICK trajectory has the following information: time (defined by start and end year), local sea level, an RCP scenario, and a BRICK dynamics model (fast or slow)

> Ruckert, K. L., Srikrishnan, V., & Keller, K. (2019). Characterizing the deep uncertainties surrounding coastal flood hazard projections: A case study for Norfolk, VA. Scientific Reports, 9(1), 1–12. https://doi.org/10/ggfsbt
"""

# ╔═╡ 9c416662-ce22-4b35-9d34-7474551811dc
struct BRICKtrajectory{I<:Integer,T<:AbstractFloat}
    years::UnitRange{I}
    lsl::Vector{<:Unitful.Length}
    rcp::T
    dynamics::AbstractString
end

# ╔═╡ b69d8c44-3271-4f15-bcd2-bab67d3807b9
md"""
Parse the original BRICK sea level simulations.
It's very quick so no need to cache.
Note that the default value for `fname` assumes the particular directory structure of this project; if you are using this for another project you will need to specify the filename explicitly or else change this argument manually.
"""

# ╔═╡ 7b384702-3938-4b71-8df5-68ddcbd8c0dd
function get_norfolk_brick(
    fname::AbstractString = datadir("raw", "brick", "BRICK.nc");
    syear::Union{Int,Missing} = missing,
    eyear::Union{Int,Missing} = missing,
)

    @assert isfile(fname) "$fname is not a valid filename"

    rcp_scenarios = NetCDF.ncread(fname, "RCP")
    dynamics = NetCDF.ncread(fname, "dynamics")
    years = Int.(NetCDF.ncread(fname, "time_proj"))
    sims = NetCDF.ncread(fname, "simulation")
    lsl_m = Array{typeof(rcp_scenarios[1]),4}(NetCDF.ncread(fname, "lsl_m"))u"m" 
	lsl = uconvert.(u"ft", lsl_m)

    nsim = length(sims)
    nrcp = length(rcp_scenarios)
    ndynam = length(dynamics)

    if ismissing(syear)
        syear = minimum(years)
    else
        @assert minimum(years) <= syear
    end
    if ismissing(eyear)
        eyear = maximum(years)
    else
        @assert eyear <= maximum(years)
    end

    idx_eyear = findfirst(years .== eyear)
    idx_syear = findfirst(years .== syear)
    trajectories = [
        BRICKtrajectory(
            syear:eyear,
            lsl[idx_syear:idx_eyear, i, j, k][:],
            rcp_scenarios[j],
            dynamics[k],
        ) for i = 1:nsim, j = 1:nrcp, k = 1:ndynam
    ][:]
    return trajectories
end;

# ╔═╡ 09815c58-3c3b-4b1d-a9a1-8a7191774019
all_trajs = get_norfolk_brick(syear=2020, eyear=2100)

# ╔═╡ c73ddb13-50d4-4e6d-ac54-c171c66b3786
begin
	function plot_lsl_models(;ub = 0.95, lb = 0.05)
		p = plot(
			legend=:topleft,
			xlabel="Year",
			ylabel="Mean Sea Level at Norfolk, VA [ft]",
			size = (500, 500),
		)
		models = ((2.6, "slow"), (4.5, "slow"), (6.0, "fast"), (8.5, "fast"))
		for (model, color) in zip(models, colors)
			rcp, dyn = model
			trajs = [
				traj for
				traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == dyn)
			]
			sims = ustrip.(u"ft", hcat([traj.lsl for traj in trajs]...))
			N = size(sims)[1]
			upper = [quantile(sims[i, :], ub) for i = 1:N]
			lower = [quantile(sims[i, :], lb) for i = 1:N]
			plot!(
				p,
				first(trajs).years,
				upper,
				fillrange = lower,
				fillalpha = 0.45,
				label = "RCP $rcp / BRICK $dyn",
				color = color,
				linewidth = 1,
			)
			plot!(p, first(trajs).years, lower, label=false, color=color, linewidth=1)
		end
		return p
	end
	p2 = plot_lsl_models()
	savefig(p2, plotsdir("msl_ranges_threemodels.pdf"))
	plot(p2)
end

# ╔═╡ 4e5b01a7-bf2c-4092-b9e2-580e46c289f5
DrWatson.wsave(
	datadir("processed", "msl_brick.jld2"),
	Dict("trajectories" => all_trajs),
)

# ╔═╡ 772a3c9f-720e-48c0-9677-2871bd4c17a0
md"Get the data from a particular year"

# ╔═╡ d0468587-2487-428c-bcb3-4e3528302147
function get_year_data(t::BRICKtrajectory, y::Int)
    idx = findfirst(t.years .== y)
    return t.lsl[idx]
end;

# ╔═╡ e969b4ce-7232-4e5b-8c28-4bbd9f62cf61
md"Get the data from a particular year"

# ╔═╡ 976185f0-2a07-4f08-b003-e5c2a58d7a97
function get_year_data(ts::Vector{<:BRICKtrajectory}, y::Int)
    idx = findfirst(first(ts).years .== y)
    return [t.lsl[idx] for t in ts]
end;

# ╔═╡ 3827f654-4a9d-4b1c-b0f7-15b7dd0e50c7
md"""
Get the NOAA sea level rise scenarios
"""

# ╔═╡ e4b61388-9e11-4220-8315-4a4d015ba71e
function get_noaa_scenarios()
	return CSV.File(datadir("raw", "noaa_sewells_point_scenarios.csv")) |> DataFrame
end;

# ╔═╡ eb8aa132-023a-4cb1-8799-87ec85621c98
noaa = get_noaa_scenarios()

# ╔═╡ 2022a380-290f-45fb-930c-e0b00d1806ee
begin
	function plot_msl_2100()

		all_trajs = get_norfolk_brick()
		noaa = get_noaa_scenarios()
		noaa_2100 = select(filter(row -> row[:year] == 2100, noaa), Not(:year))

		p = plot(
			xlabel = "MSL at Norfolk, VA in 2100 [ft]",
			ylabel = "Cumulative Probability",
			size = (500, 500),
		)

		all_rcp = unique([traj.rcp for traj in all_trajs])
		all_dyn = unique([traj.dynamics for traj in all_trajs])
		all_msl_2100 = [
			ustrip(u"ft", get_year_data(traj, 2100)) for traj in all_trajs
		]

		for (rcp, color) in zip(all_rcp, colors[length(all_rcp):end])
			for (dyn, lsty) in zip(all_dyn, [:solid, :dash])
				trajs = [
					traj for
					traj in all_trajs if (traj.rcp == rcp) & (traj.dynamics == dyn)
				]
				sims = ustrip.(u"ft", [get_year_data(traj, 2100) for traj in trajs])
				label = "RCP $rcp / $dyn"
				empirical = ecdf(sims)
				plot!(
					p,
					empirical,
					linestyle=lsty,
					linewidth=2,
					label=label,
					color=color,
				)
			end
		end
		
		for scenario in names(noaa_2100)
			y = noaa_2100[1, scenario]
			Plots.vline!(
				p,
				[y],
				color=:gray,
				label=false,
				linestyle=:dot,
				linewidth=2,
			)
			annotate!(
				p,
				y + 0.25,
				0.05,
				text("NOAA $scenario", :right, 6, color = :gray, rotation=270),
			)
		end
		
		plot!(p, legend=:right)

		return p
	end
	p1 = plot_msl_2100()
	savefig(p1, plotsdir("msl_multiple_pdfs.pdf"))
	plot(p1)
end

# ╔═╡ eece84cd-3f44-4d8b-851e-933d22151c19
begin
	function plot_all_trajs()
		all_trajs = get_norfolk_brick(syear=2000, eyear=2100)
		noaa = get_noaa_scenarios()
		p = plot(
			xlabel = "Year",
			ylabel = "Mean Sea Level at Norfolk, VA [ft]",
			legend = :topleft,
			size = (500, 500),
		)
		for traj in sample(all_trajs, 2_500)
			plot!(
				p,
				traj.years,
				ustrip.(u"ft", traj.lsl),
				linewidth=0.25,
				alpha=0.25,
				color=:gray,
				label=false,
			)
		end
		for (scenario, color) in zip(names(select(noaa, Not(:year))), colors)
			plot!(
				noaa[!, "year"],
				noaa[!, scenario],
				color=color,
				label="NOAA $scenario",
				linewidth=1,
			)
		end
		return p
	end
	p3 = plot_all_trajs()
	savefig(p3, plotsdir("msl_all_trajs.png"))
	p3
end

# ╔═╡ 84182c9e-1b97-4e3e-b3dc-b1c8c50d899c
begin
	function make_boxplots()

		all_trajs = get_norfolk_brick()
		noaa = get_noaa_scenarios()
		
		rows = [
			DataFrame(
				"scenario" => "RCP $(traj.rcp), $(traj.dynamics) dynamics",
				"msl_2100" => ustrip(u"ft", get_year_data(traj, 2100))
			)
			for traj in all_trajs
		]
		lsl_2100 = vcat(rows...)

		p = plot(
			ylabel = "MSL at Norfolk, VA in 2100 [ft]",
			bottom_margin = 10mm,
			left_margin = 5mm,
			size = (700, 300),
			xtickfontsize = 7,
			ylabelfontsize = 9,
		)
		@df lsl_2100 violin!(p, :scenario, :msl_2100, xrotation=30, legend=false)
		
		return p
	end
	p4 = make_boxplots()
	savefig(p4, plotsdir("msl_boxplots.pdf"))
	plot(p4)
end

# ╔═╡ Cell order:
# ╟─49cd8998-378d-11ec-3ceb-811bf9c58697
# ╠═6eaeeaaf-5fe8-44b2-8957-b47918f26b17
# ╟─ab567401-943f-4877-8572-cfa8ffb531d4
# ╟─c259b7d3-0873-4502-9d4b-a5d18c5adcd8
# ╠═09815c58-3c3b-4b1d-a9a1-8a7191774019
# ╠═eb8aa132-023a-4cb1-8799-87ec85621c98
# ╟─5048e15a-d382-4a3c-9aa6-d0d34cb2b297
# ╠═a4a6d097-0d5c-48c7-b8bc-34baa3c910cc
# ╟─f33c8247-c110-412d-a7d4-4b8aa8f89f34
# ╟─2022a380-290f-45fb-930c-e0b00d1806ee
# ╟─37d98257-5de6-456a-a542-36d0dc1c7a96
# ╟─f251e057-4409-4c8b-a2d7-4603c8af6750
# ╟─c73ddb13-50d4-4e6d-ac54-c171c66b3786
# ╠═5f912aeb-c9e4-485a-897f-a71e84ec6038
# ╟─eece84cd-3f44-4d8b-851e-933d22151c19
# ╟─f815e33f-97fa-4a1b-b8fc-77241d81d172
# ╠═84182c9e-1b97-4e3e-b3dc-b1c8c50d899c
# ╟─48b37aa3-b756-4628-be82-5289ed202aa2
# ╠═4e5b01a7-bf2c-4092-b9e2-580e46c289f5
# ╟─15ad5394-ee85-4aff-851d-e9e9a58ff1a0
# ╟─3a5a60eb-3493-4359-924e-ad19ddd28b63
# ╠═9c416662-ce22-4b35-9d34-7474551811dc
# ╟─b69d8c44-3271-4f15-bcd2-bab67d3807b9
# ╠═7b384702-3938-4b71-8df5-68ddcbd8c0dd
# ╟─772a3c9f-720e-48c0-9677-2871bd4c17a0
# ╠═d0468587-2487-428c-bcb3-4e3528302147
# ╟─e969b4ce-7232-4e5b-8c28-4bbd9f62cf61
# ╠═976185f0-2a07-4f08-b003-e5c2a58d7a97
# ╠═3827f654-4a9d-4b1c-b0f7-15b7dd0e50c7
# ╠═e4b61388-9e11-4220-8315-4a4d015ba71e
