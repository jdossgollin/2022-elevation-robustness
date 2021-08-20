### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ d32837b2-e05e-11eb-0fb7-1fc68c7f97bf
begin
	using ColorSchemes
	using Distributions
	using DrWatson
	using KernelDensity
	using PlutoUI
	using Plots
	using StatsBase
	using StatsPlots
	using Unitful
	using Plots: mm

    using NorfolkBRICK

    PlutoUI.TableOfContents()
end

# ╔═╡ cbc1f922-4e33-4bf9-bb2c-ee93e7086c68
md"""
# Probabilistic MSL Analysis
"""

# ╔═╡ 50c162f8-b90a-43a6-8c9e-0074329dbbe6
colors = ColorSchemes.okabe_ito;

# ╔═╡ c981f860-7dbe-4313-8d6f-7268b4249ec7
md"""
We combine melding and rejection sampling
"""

# ╔═╡ a52e2861-1717-424d-b8f7-1d43f6c8efac
trajectories = get_norfolk_brick()

# ╔═╡ 7af03bed-6341-48bb-b04d-ac9f584af58b
get_lsl_2100(trajs) = ustrip.(u"ft", get_year_data(trajs, 2100));

# ╔═╡ f2c4ee04-04a8-4eaa-bfae-082b850c685d
function normalize(w::Vector{<:Real})
	return StatsBase.pweights(w ./ sum(w))
end

# ╔═╡ 315a762a-fce2-444b-bf2d-395d3b736457
function rejection_sampling(trajs, target)
	y = get_lsl_2100(trajs)
	target_trunc = truncated(target, quantile(y, 0.0025), quantile(y, 0.9975))
	sampling = kde(y)
	w = pdf(target_trunc, y) ./ pdf(sampling, y)
	return w
end

# ╔═╡ 9a83aa2b-f2cd-4c73-a2a5-20fbed0278e7
function get_unit_weights(trajs::Vector{<:BRICKtrajectory})
	N = length(trajs)
	w = ones(N) ./ N
	return normalize(w)
end

# ╔═╡ 6a361969-3851-44b6-b2fe-53f295cc2dde
function get_subjective_weights(trajs::Vector{<:BRICKtrajectory})
	# https://coast.noaa.gov/slr/#/layer/sce/0/-8491236.598143658/4473859.265531366/9/satellite/22/0.8/2050/inter/midAccretion
	target = Normal(4, 1)
	w_pdf = rejection_sampling(trajs, target)
	
	w_rcp = Dict(2.6 => 1, 4.5 => 3, 6.0 => 3, 8.5 => 1)
	w_dyn = Dict("fast" => 0.5, "slow" => 0.5)
	w_physics = [w_rcp[traj.rcp] * w_dyn[traj.dynamics] for traj in trajs]
	
	mixing = 0.99
	w_melded =  @. (w_pdf ^ mixing) * (w_physics ^ (1 - mixing))
	return normalize(w_melded)
end

# ╔═╡ 692cea94-9b2f-432b-b210-77942b0c2084
function get_uncertain_weights(trajs::Vector{<:BRICKtrajectory})
	target = Beta(2, 5)
	y = get_lsl_2100(trajs)
	ub = maximum(y)
	lb = minimum(y)
	rescaled = @. (y - lb) / (ub - lb)
	sampling = kde(rescaled)
	w = pdf(target, rescaled) ./ pdf(sampling, rescaled)
	return normalize(w)
end

# ╔═╡ e32d46b3-a358-464f-a79c-fe5ea3d45de6
function get_confident_weights(trajs::Vector{<:BRICKtrajectory})
	w_rcp = Dict(2.6 => 1, 4.5 => 5, 6.0 => 5, 8.5 => 1)
	w_dyn = Dict("fast" => 0.5, "slow" => 0.5)
	w = [w_rcp[traj.rcp] * w_dyn[traj.dynamics] for traj in trajs]
	return normalize(w)
end

# ╔═╡ b8780bf2-178c-4f9f-9084-8dcd94bb8b66
function get_naive_weights(trajs::Vector{<:BRICKtrajectory})
	target = Uniform(1, 5)
	w = rejection_sampling(trajs, target)
	return normalize(w)
end

# ╔═╡ fbdae066-51a9-4e49-a7d8-5c3676cff61d
ww = Dict(
	"1: Do Not Weight" => get_unit_weights(trajectories),
	"2: My Best Guess" => get_subjective_weights(trajectories),
	"3: High Certainty" => get_confident_weights(trajectories),
	"4: Low Certainty" => get_uncertain_weights(trajectories),
	"5: Naive Uniform" => get_naive_weights(trajectories),
);

# ╔═╡ e21c84ce-6251-4c3b-813f-e1f01b73dab6
begin
	p3 = plot(
		xlabel = "Mean Sea Level in 2100 [ft]",
		ylabel = "Probability Density",
		yticks = [],
		left_margin = 5mm,
	)
	for (model, c, lst) in zip(keys(sort(ww)), colors, 2:10)
		density!(
			p3,
			get_lsl_2100(trajectories),
			weights = ww[model],
			label = model,
			normalize=true,
			color = c,
			linewidth = 2,
			linestyle = lst,
		)
	end
	p3
end

# ╔═╡ Cell order:
# ╠═d32837b2-e05e-11eb-0fb7-1fc68c7f97bf
# ╠═cbc1f922-4e33-4bf9-bb2c-ee93e7086c68
# ╠═50c162f8-b90a-43a6-8c9e-0074329dbbe6
# ╠═c981f860-7dbe-4313-8d6f-7268b4249ec7
# ╠═a52e2861-1717-424d-b8f7-1d43f6c8efac
# ╠═7af03bed-6341-48bb-b04d-ac9f584af58b
# ╠═f2c4ee04-04a8-4eaa-bfae-082b850c685d
# ╠═315a762a-fce2-444b-bf2d-395d3b736457
# ╠═9a83aa2b-f2cd-4c73-a2a5-20fbed0278e7
# ╠═6a361969-3851-44b6-b2fe-53f295cc2dde
# ╠═692cea94-9b2f-432b-b210-77942b0c2084
# ╠═e32d46b3-a358-464f-a79c-fe5ea3d45de6
# ╠═b8780bf2-178c-4f9f-9084-8dcd94bb8b66
# ╠═fbdae066-51a9-4e49-a7d8-5c3676cff61d
# ╠═e21c84ce-6251-4c3b-813f-e1f01b73dab6
