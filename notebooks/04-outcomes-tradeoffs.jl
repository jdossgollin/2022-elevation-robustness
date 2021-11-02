### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 62956480-37bc-11ec-1901-13f49754cd1d
begin
	using Pkg
	Pkg.activate("..")
	using ColorSchemes
	using CSV
	using DataFrames
	using Distributions
	using Downloads
	using DrWatson
	using KernelDensity
	using PlutoUI
	using StaticArrays
	using StatsBase
	using StatsPlots
	using Unitful
	using Turing
	TableOfContents()
end

# ╔═╡ 3c2d6a02-1f68-48ef-9c9c-3b29a653b2c5
md"""
# Notebook 04: Outcomes and Tradeoffs

> Doss-Gollin and Keller, 2021, In Prep. **Do not circulate!**

*If you have not already run Notebooks 01, 02, and 03 stop and do that!*
"""

# ╔═╡ c57017ea-d32d-4fb0-8cf7-88ffc45edc93
colors = ColorSchemes.okabe_ito;

# ╔═╡ 24a51c14-9046-42d0-8999-2b0c5e6ee519
md"""
## Load data

We start by loading the scenarios (𝐬; `\bfs`), decisions (𝐱; `\bfx`), and outcomes (y).
"""

# ╔═╡ d92e3059-9bc5-4993-a38d-6dfad3854ae5
begin
	cachefile = datadir("processed", "sows_outcomes.jld2")
	u = DrWatson.load(cachefile, "u")
	𝐬 = DrWatson.load(cachefile, "𝐬")
	𝐱 = DrWatson.load(cachefile, "𝐱")
end;

# ╔═╡ 4544a2aa-0926-45c9-9082-5f706f7b9273
md"""
## Bayesian weighting

>TODO add some explanation
"""

# ╔═╡ 37b2cda6-7d57-404a-a235-fb4ec04f26c4
begin
	lsl_2100 = ustrip.(u"ft", map(s -> s.lsl[findfirst(s.years .== 2100)], 𝐬))
	lb = minimum(lsl_2100)#quantile(lsl_2100, 0.0025)
	ub = maximum(lsl_2100)#quantile(lsl_2100, 0.9975)
	(lb, ub)
end

# ╔═╡ 9266995f-9a66-4bf1-acbf-4ee515cd2900
begin
	priors = [
		LogNormal(0.4, 0.4),
		LogNormal(1, 0.35),
		LogNormal(0.9, 0.5),
	]
	prior_names = ["Slow SLR", "Fast SLR", "Uncertain SLR"]
end;

# ╔═╡ 7c58c2a0-87df-4a90-8d41-1bca62c732db
md"""
## Prior over model structures and parameters
"""

# ╔═╡ 88853a76-196b-426a-98c4-6c4d038b5dfa
begin
	function plot_priors()
		p = plot(
			xlabel = "LSL in 2100",
			ylabel = "Probability Density",
			xlims = (0, 8),
			size = (500, 500),
		)
		for (i, prior) in enumerate(priors)
			plot!(p, prior, label=prior_names[i], c=colors[i])
		end
		density!(p, lsl_2100, label="Uniform Weights", c=colors[length(priors)+1])
		return p
	end
	p_prior = plot_priors()
	savefig(p_prior, plotsdir("lsl_priors.pdf"))
	p_prior
end

# ╔═╡ 4bfe96ce-0a5b-4496-b55e-ba03b95f39dc
md"""
### Diagnostics

It's helpful to see what the weights are that we're assigning to different scenarios.
We don't want any single scenario to have *too* much leverage...
"""

# ╔═╡ a641f9cc-dd1e-409f-8ad6-1afe78f3012a
md"""
## Tradeoff plots

Now we get to the fun stuff!!!
"""

# ╔═╡ 6e81c5b2-ede4-4197-aa29-a41a3dfcae34
md"""
## Helper functions
"""

# ╔═╡ 13edb614-2d68-4c4e-beb8-22c916b7cccd
md"Normalize weights so they sum to 1"

# ╔═╡ 2ae53ed4-1671-4fb6-a13e-9505d02d262a
function normalize(w::Vector{<:Real})
	return StatsBase.pweights(w ./ sum(w))
end;

# ╔═╡ 60d96657-91db-4856-8677-aa388f6d59af
scenario_weights = [
	(
		name="RCP $rcp / BRICK $dynamics",
		w=normalize([
			float(Int((s.rcp == rcp) & (s.dynamics == dynamics)))
			for s in 𝐬
		]),
	)
	for (rcp, dynamics) in zip([2.6, 4.5, 6.0, 8.5], ["slow", "slow", "fast", "fast"])
];

# ╔═╡ 280c18b6-82bf-4eb2-b77d-6fc53557e84c
md"Calculate the 95th percentile minus 5th percentile, with weights"

# ╔═╡ 2f914d70-a6eb-4f2c-8b08-403ecef72a5a
function qrange(x, w; ub=0.9, lb=0.1)
	upper = quantile(x, w, ub)
	lower = quantile(x, w, lb)
	return upper - lower
end;

# ╔═╡ 423b133f-8b46-444d-8db3-da1df584e35a
md"""
For a given weight over each scenario, calculate metrics.
These are indexed [decision, metric].
"""

# ╔═╡ 4ad926d8-484e-43d9-93a4-ff0709ff1361
function calc_metrics(w)
	
	I, J, K = size(u)
	@assert length(w) == I

	M = 5 # number of metrics to calculate
	metrics = zeros(Float64, J, M)
	metric_names = repeat([""], M)
	
	# first metric: height increase in feet
	metrics[:, 1] = ustrip.(u"ft", 𝐱)
	metric_names[1] = "Height Increase [ft]"

	# second metric: construction cost in USD
	metrics[:, 2] = u[1, :, 1]
	metric_names[2] = "Construction Cost [USD]"

	# third metric: expected NPV damages in USD
	metrics[:, 3] = [mean(u[:, j, 2], w) for j in 1:J]
	metric_names[3] = "Expected NPV Lifetime Flood Damages [USD]"

	# fourth metric: expected lifetime costs in USD
	metrics[:, 4] = [mean(u[:, j, 3], w) for j in 1:J]
	metric_names[4] = "Expected NPV Lifetime Total Costs [USD]"

	# fifth metric: 90% CI lifetime costs
	metrics[:, 5] = [qrange(u[:, j, 3], w) for j in 1:J]
	metric_names[5] = "80% CI Width of NPV Lifetime Total Costs [USD]"
	
	return metrics, metric_names
end;

# ╔═╡ 26525f5f-63bb-4aec-9078-bc7c509876ae
md"""
Create a generic function to plot 2D tradeoffs
"""

# ╔═╡ 2119c4c3-4ee4-479a-aed5-7477b4980507
function plot_tradeoff(weights, kx, ky; kwargs...)
	_, mnames = calc_metrics(first(weights).w)
	p = plot(
		xlabel = mnames[kx],
		ylabel = mnames[ky];
		kwargs...
	)
	for (color, weight) in zip(colors, weights)
		metrics, _ = calc_metrics(weight.w)
		m̂x = metrics[:, kx]
		m̂y = metrics[:, ky]
		plot!(p, m̂x, m̂y, color=color, label=weight.name, marker="-.")
	end
	return p
end;

# ╔═╡ 95add889-a398-4bae-a3cb-37ca7c620644
begin
	pt3 = plot_tradeoff(scenario_weights, 1, 4; size=(500, 500), legend=:bottomright)
	savefig(pt3, plotsdir("tradeoffs_height_totalcost_byscenario.pdf"))
	pt3
end

# ╔═╡ 37f159ec-1cd2-4463-8862-e39bae53e9cb
begin
	pt4 = plot_tradeoff(scenario_weights, 1, 3; size=(500, 500), legend=:topright)
	savefig(pt4, plotsdir("tradeoffs_height_damages_byscenario.pdf"))
	pt4
end

# ╔═╡ 1fd21800-92cb-4b8b-bd4b-3c4afa73b43a
begin
	pt6 = plot_tradeoff(scenario_weights, 1, 5; size=(500, 500), legend=:topright)
	savefig(pt6, plotsdir("tradeoffs_height_iqr_byscenario.pdf"))
	pt6
end

# ╔═╡ 46a288ee-3c29-497d-964e-7aeb534ab886
md"""
Function to get some weights
>TODO add doc
"""

# ╔═╡ 675e0f3f-988a-4541-8e3e-c3443365ac6b
function get_weights(svec, d)
	target_dist = truncated(d, lb, ub)
	y = ustrip.(u"ft", map(s -> s.lsl[findfirst(s.years .== 2100)], svec))
	sampling_dist = kde(y)
	p = pdf(target_dist, y) ./ pdf(sampling_dist, y)
	w = normalize(p)
	return w
end;

# ╔═╡ cda94344-47c8-437a-a14c-cc475d66e096
begin
	weights = [get_weights(𝐬, d) for d in priors]
	push!(weights, normalize(repeat([1], length(𝐬))))
	push!(prior_names, "Uniform Weights")
end;

# ╔═╡ 4b498f43-0617-45f8-8ae3-f93fde1cfa83
prior_weights = [(name=n, w=w) for (n, w) in zip(prior_names, weights)];

# ╔═╡ 31620a55-ce09-4de7-a5b5-b2f6ef850abd
begin
	pt1 = plot_tradeoff(prior_weights, 1, 4; size=(500, 500), legend=:bottomright)
	savefig(pt1, plotsdir("tradeoffs_height_totalcost_byprior.pdf"))
	pt1
end

# ╔═╡ 4147df61-2b35-4583-ba5e-c55c6f90d5cb
begin
	pt2 = plot_tradeoff(prior_weights, 1, 3; size=(500, 500), legend=:topright)
	savefig(pt2, plotsdir("tradeoffs_height_damages_byprior.pdf"))
	pt2
end

# ╔═╡ 5bead18e-ee79-46ad-9cbe-b234791561b8
begin
	pt5 = plot_tradeoff(prior_weights, 1, 5; size=(500, 500), legend=:topright)
	savefig(pt5, plotsdir("tradeoffs_height_iqr_byprior.pdf"))
	pt5
end

# ╔═╡ 6bfb6f34-4499-4db3-8038-bfab3e8d4756
begin
	function plot_scenario_weights()
		yticks = [1/1000, 1/100, 1/10, 1, 10, 100]
		p = plot(
			xlabel="MSL in 2100",
			ylabel="Relative Scenario Weight",
			yscale=:log10,
			legend=:top,
			yticks = (yticks, string.(yticks)),
			dpi = 200,
		)
		for (i, w) in enumerate(weights)
			scatter!(
				lsl_2100,
				w ./ mean(w),
				color=colors[i],
				label=prior_names[i],
				markerstrokewidth=0,
				markersize=2,
			)
		end
		return p
	end
	p_scenario_weights = plot_scenario_weights()
	savefig(p_scenario_weights, plotsdir("scenario_weights.png"))
	p_scenario_weights
end

# ╔═╡ Cell order:
# ╟─3c2d6a02-1f68-48ef-9c9c-3b29a653b2c5
# ╠═62956480-37bc-11ec-1901-13f49754cd1d
# ╠═c57017ea-d32d-4fb0-8cf7-88ffc45edc93
# ╟─24a51c14-9046-42d0-8999-2b0c5e6ee519
# ╠═d92e3059-9bc5-4993-a38d-6dfad3854ae5
# ╟─4544a2aa-0926-45c9-9082-5f706f7b9273
# ╠═37b2cda6-7d57-404a-a235-fb4ec04f26c4
# ╠═9266995f-9a66-4bf1-acbf-4ee515cd2900
# ╟─7c58c2a0-87df-4a90-8d41-1bca62c732db
# ╟─88853a76-196b-426a-98c4-6c4d038b5dfa
# ╠═cda94344-47c8-437a-a14c-cc475d66e096
# ╠═4b498f43-0617-45f8-8ae3-f93fde1cfa83
# ╠═60d96657-91db-4856-8677-aa388f6d59af
# ╟─4bfe96ce-0a5b-4496-b55e-ba03b95f39dc
# ╟─6bfb6f34-4499-4db3-8038-bfab3e8d4756
# ╟─a641f9cc-dd1e-409f-8ad6-1afe78f3012a
# ╠═31620a55-ce09-4de7-a5b5-b2f6ef850abd
# ╠═4147df61-2b35-4583-ba5e-c55c6f90d5cb
# ╟─5bead18e-ee79-46ad-9cbe-b234791561b8
# ╟─95add889-a398-4bae-a3cb-37ca7c620644
# ╟─37f159ec-1cd2-4463-8862-e39bae53e9cb
# ╟─1fd21800-92cb-4b8b-bd4b-3c4afa73b43a
# ╟─6e81c5b2-ede4-4197-aa29-a41a3dfcae34
# ╟─13edb614-2d68-4c4e-beb8-22c916b7cccd
# ╠═2ae53ed4-1671-4fb6-a13e-9505d02d262a
# ╟─280c18b6-82bf-4eb2-b77d-6fc53557e84c
# ╠═2f914d70-a6eb-4f2c-8b08-403ecef72a5a
# ╠═423b133f-8b46-444d-8db3-da1df584e35a
# ╠═4ad926d8-484e-43d9-93a4-ff0709ff1361
# ╟─26525f5f-63bb-4aec-9078-bc7c509876ae
# ╠═2119c4c3-4ee4-479a-aed5-7477b4980507
# ╟─46a288ee-3c29-497d-964e-7aeb534ab886
# ╠═675e0f3f-988a-4541-8e3e-c3443365ac6b
