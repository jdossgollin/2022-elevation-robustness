### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ a5a31f93-c222-4478-a077-f6cfbb660dc2
using Pkg; Pkg.add("Optim")

# ╔═╡ ca64c816-a7bc-11eb-38c7-97aa557f6992
using DrWatson

# ╔═╡ c4d1b067-94f9-4f0d-8bcc-5b2a61d439e5
begin
	using PlutoUI;
	PlutoUI.TableOfContents();
end

# ╔═╡ 80d8bfb5-984c-48a3-b682-b7c0ad6ac13c
begin
    using ColorSchemes
    using Distributions
	using Optim
	using Plots
	using StatsBase
	using Turing
    using Unitful
    using UnitfulRecipes

    using DecisionProblem
    using HouseElevation
    using NorfolkFloods

    colors = ColorSchemes.tab10
end;

# ╔═╡ b6e06c8a-ef9c-4ec4-88ba-bd3c9f36b237
@quickactivate "2021-elevation-robustness"

# ╔═╡ 109d98db-71c3-425e-b69b-2ee54e244ce8
md"""
# Analysis
"""

# ╔═╡ 9e6daa36-6343-49f4-b2b6-1df2755617e9
md"Load the SOWs from file"

# ╔═╡ fcbac66b-2e4e-4a74-8571-b4d2efd65ce4
sows = DrWatson.load(datadir("processed", "SOWs.jld2"), "sows");

# ╔═╡ e5a57393-b57d-4a11-9a72-2f992501d780
md"""
## Calculate BFE

Zarekarizi et al (2020) defines the house elevation with respect to the BFE, which we define as the stationary 100 year flood
"""

# ╔═╡ f2bf864b-937e-4d8f-b84e-4414d9cb5474
obs = NorfolkFloods.get_norfolk_annual();

# ╔═╡ a494aace-99c7-4535-853a-6fb8ed059709
@model function GEV(y)
    μ ~ Turing.Flat()
    σ ~ Turing.FlatPos(0.025)
    ξ ~ Turing.Flat()
    y .~ Distributions.GeneralizedExtremeValue(μ, σ, ξ)
end;

# ╔═╡ 7cfa3c22-cdcd-48b6-b96d-bd3500076b78
map_estimate = optimize(GEV(ustrip.(u"ft", obs.surge)), MAP(), NelderMead())

# ╔═╡ e4c3434e-4114-47bd-8e36-c2696f1c3aca
map_gev = GeneralizedExtremeValue(3.62221, 0.503049, 0.203107); # manual = consistent

# ╔═╡ bc17ce92-9e10-47e0-ac63-0e49606e6083
surge_q99 = quantile(map_gev, 0.99)u"ft"

# ╔═╡ be8d5b6d-d2c6-460c-89bb-cf8d7010648d
bfe = uconvert(u"ft", surge_q99 + maximum(obs.msl))

# ╔═╡ a3362135-b5e1-4b9b-a610-2409885d4f3a
md"This tells us that in the year 2020, we would expect that a flood exceed the BFE ($(bfe)) with probability 1/100. This allows us to define our house:"

# ╔═╡ 4ef02fd2-430b-438c-85fd-c75eaa284476
house = HouseStructure(1000u"ft^2", 250_000, bfe - 1u"ft")

# ╔═╡ 4ecf8863-8b05-46e2-92e9-8a9168427b3c
md"""
## Evaluate SOWs

First, we evaluate the SOWs
"""

# ╔═╡ 30136674-10a3-4f7a-8ddb-d0088d1376fd
Δh = (0:7:14)u"ft";

# ╔═╡ 9433a43e-67e2-45aa-9cbc-8aeab464f7f7
γ = 1.0 - 0.03

# ╔═╡ 20113e6e-fa97-4cff-a812-da5001c9e22e
outcomes = [evaluate(house, Lever(hh), sows; γ = γ) for hh in Δh]

# ╔═╡ 31c65a1d-c6ba-4b58-9f57-57caafb43f18
md"""
## Combine

We can combine by weighting
"""

# ╔═╡ d7ecebb8-07c6-4694-a4c4-389e8bfb6815
function combine(metrics, w)
	avgs = map(x -> mean(x, w), eachcol(metrics))
    return DataFrame([name => avg for (name, avg) in zip(names(metrics), avgs)]...)
end

# ╔═╡ 1bdbc101-8aa9-4190-8009-58ee8cfe8b3b
md"### Uniform Weights"

# ╔═╡ 5f65ff9f-caff-4adb-af9d-bbe5f8a65764
w_unif = StatsBase.uweights(length(sows.sows));

# ╔═╡ c4e2e4f8-879f-4266-8f6f-eaf897c4f7b7
avgs_unif = vcat([combine(oc, w_unif) for oc in outcomes]...)

# ╔═╡ 54f787e6-d748-47f2-a37c-6b354385f445
function plot_results(results, hh)
	p = plot(
		hh,
		results[:construction_cost],
		label="Construction Cost",
		xlabel="House Elevation in 2020",
		ylabel="NPV Costs (USD2020)",
		legend=(0.7, 0.3),
		w=2,
	)
	plot!(p, hh, results[:npv_damages], label="NPV Damages", w=2)
	plot!(p, hh, results[:total_cost], label="Total Cost", w=3)
	best_hgt = argmin(results[:total_cost])
	scatter!(p, [hh[best_hgt]], [results[:total_cost][best_hgt]], label="Optimal Height", m=7)
	return p
end

# ╔═╡ 7a08cbc2-0e1d-4932-a924-9f98e9ce9fa4
plot(plot_results(avgs_unif, Δh), title="All Scenarios Weighted Equal; House $(house.h)")

# ╔═╡ 0804d98d-4310-4adb-ac8f-7dd729682065
md"""
## Specific Scenario
"""

# ╔═╡ eadaa7ce-75c6-4f3e-936b-ebea99f2dc7b
sow = first(sows.sows)

# ╔═╡ ada38d3a-8adc-4b01-a4d4-5d0654c10308
msl_2100 = 

# ╔═╡ 4de0a20d-0acf-4f75-9583-3ce08a90100d


# ╔═╡ Cell order:
# ╠═ca64c816-a7bc-11eb-38c7-97aa557f6992
# ╠═b6e06c8a-ef9c-4ec4-88ba-bd3c9f36b237
# ╠═c4d1b067-94f9-4f0d-8bcc-5b2a61d439e5
# ╠═109d98db-71c3-425e-b69b-2ee54e244ce8
# ╠═a5a31f93-c222-4478-a077-f6cfbb660dc2
# ╠═80d8bfb5-984c-48a3-b682-b7c0ad6ac13c
# ╠═9e6daa36-6343-49f4-b2b6-1df2755617e9
# ╠═fcbac66b-2e4e-4a74-8571-b4d2efd65ce4
# ╠═e5a57393-b57d-4a11-9a72-2f992501d780
# ╠═f2bf864b-937e-4d8f-b84e-4414d9cb5474
# ╠═a494aace-99c7-4535-853a-6fb8ed059709
# ╠═7cfa3c22-cdcd-48b6-b96d-bd3500076b78
# ╠═e4c3434e-4114-47bd-8e36-c2696f1c3aca
# ╠═bc17ce92-9e10-47e0-ac63-0e49606e6083
# ╠═be8d5b6d-d2c6-460c-89bb-cf8d7010648d
# ╟─a3362135-b5e1-4b9b-a610-2409885d4f3a
# ╠═4ef02fd2-430b-438c-85fd-c75eaa284476
# ╠═4ecf8863-8b05-46e2-92e9-8a9168427b3c
# ╠═30136674-10a3-4f7a-8ddb-d0088d1376fd
# ╠═9433a43e-67e2-45aa-9cbc-8aeab464f7f7
# ╠═20113e6e-fa97-4cff-a812-da5001c9e22e
# ╠═31c65a1d-c6ba-4b58-9f57-57caafb43f18
# ╠═d7ecebb8-07c6-4694-a4c4-389e8bfb6815
# ╠═1bdbc101-8aa9-4190-8009-58ee8cfe8b3b
# ╠═5f65ff9f-caff-4adb-af9d-bbe5f8a65764
# ╠═c4e2e4f8-879f-4266-8f6f-eaf897c4f7b7
# ╠═54f787e6-d748-47f2-a37c-6b354385f445
# ╠═7a08cbc2-0e1d-4932-a924-9f98e9ce9fa4
# ╠═0804d98d-4310-4adb-ac8f-7dd729682065
# ╠═eadaa7ce-75c6-4f3e-936b-ebea99f2dc7b
# ╠═ada38d3a-8adc-4b01-a4d4-5d0654c10308
# ╠═4de0a20d-0acf-4f75-9583-3ce08a90100d
