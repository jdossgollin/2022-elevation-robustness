### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 6eaeeaaf-5fe8-44b2-8957-b47918f26b17
begin
	# we start by loading in the packages we use
	using Pkg
	Pkg.activate("..")
	using ColorSchemes
	using DataFrames
	using Distributions
	using DrWatson
	using PlutoUI
	using Plots
	using Unitful
	using Turing
	TableOfContents()
end

# ╔═╡ 49cd8998-378d-11ec-3ceb-811bf9c58697
md"""
# Notebook 03: SOWs and OUtcomes

> Doss-Gollin and Keller, 2021, In Prep. **Do not circulate!**

In this notebook we'll explore simulations of mean sea level.

*If you have not already run Notebooks 01 and 02, stop and do that!*
"""

# ╔═╡ ab567401-943f-4877-8572-cfa8ffb531d4
md"""
## Read Data

We start by loading in the previously defined data
"""

# ╔═╡ 86eb1175-f48e-456a-a4d4-1ec5d0c0a55e
posterior_df = DrWatson.load(datadir("processed", "surge_posterior.jld2"), "posterior") |> DataFrame

# ╔═╡ 04640130-791c-4aac-8af2-450fbcbf7bdb
lsl_trajs = DrWatson.load(datadir("processed", "msl_brick.jld2"), "trajectories")

# ╔═╡ 4450ffd7-679b-49f0-8128-8f2ff2beb906
md"""
## Create SOWs

Now we need to create a database of SOWs
"""

# ╔═╡ 5d8adb71-5e3a-41dd-aaad-c76cf3d43c61
n_SOWs = 100_000

# ╔═╡ d70e47ab-33d8-4e63-ac5d-45c126656d5d
SOWs

# ╔═╡ f180e72b-d2c5-46d7-a2ef-84c012e59cf4
md"""
## Outcomes
"""

# ╔═╡ fedc0001-84b8-4b2d-856d-34ac7fa99d41
md"""
### Decisions
"""

# ╔═╡ c07c39ea-da4b-4a2e-bfa2-99d4a078c513
𝐱 = (0:0.5:16)u"ft"

# ╔═╡ 26eb91da-a71b-4eb9-9a8c-a9a336b34446
md"""
### System Model
"""

# ╔═╡ e5456332-766a-48e3-82c9-34224b7589a2
md"""
## Helper Functions
"""

# ╔═╡ 16bb254d-3773-4820-b4dd-fe114f7016a4
struct SOW{I<:Integer,T<:AbstractFloat}
    years::UnitRange{I}
    lsl::Vector{<:Unitful.Length}
	gev::typeof(GeneralizedExtremeValue(1, 2, 0.))
    rcp::T
    dynamics::AbstractString
end

# ╔═╡ 89f27639-7b5f-408c-80e7-7cb16cb688c4
begin
	𝐬 = []
	for _ in 1:n_SOWs
		row = rand(1:nrow(posterior_df))
		surge = GeneralizedExtremeValue(
			posterior_df[row, :μ],
			posterior_df[row, :σ],
			posterior_df[row, :ξ],
		)
		lsl = rand(lsl_trajs)
		sow = SOW(
			lsl.years,
			lsl.lsl,
			surge,
			lsl.rcp,
			lsl.dynamics,
		)
		push!(𝐬, sow)
	end
end

# ╔═╡ Cell order:
# ╠═49cd8998-378d-11ec-3ceb-811bf9c58697
# ╠═6eaeeaaf-5fe8-44b2-8957-b47918f26b17
# ╠═ab567401-943f-4877-8572-cfa8ffb531d4
# ╠═86eb1175-f48e-456a-a4d4-1ec5d0c0a55e
# ╠═04640130-791c-4aac-8af2-450fbcbf7bdb
# ╠═4450ffd7-679b-49f0-8128-8f2ff2beb906
# ╠═5d8adb71-5e3a-41dd-aaad-c76cf3d43c61
# ╠═89f27639-7b5f-408c-80e7-7cb16cb688c4
# ╠═d70e47ab-33d8-4e63-ac5d-45c126656d5d
# ╠═f180e72b-d2c5-46d7-a2ef-84c012e59cf4
# ╠═fedc0001-84b8-4b2d-856d-34ac7fa99d41
# ╠═c07c39ea-da4b-4a2e-bfa2-99d4a078c513
# ╠═26eb91da-a71b-4eb9-9a8c-a9a336b34446
# ╠═e5456332-766a-48e3-82c9-34224b7589a2
# ╠═16bb254d-3773-4820-b4dd-fe114f7016a4
