### A Pluto.jl notebook ###
# v0.17.0

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
	using Distributions
	using DrWatson
	using Interpolations
	using MCMCChains: Chains
	using Unitful
	using UnitfulRecipes
	using PlutoUI
	using Plots
	using Tullio
	TableOfContents()
end

# ╔═╡ 49cd8998-378d-11ec-3ceb-811bf9c58697
md"""
# Notebook 03: SOWs and Outcomes

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
posterior_df = read(datadir("processed", "surge_posterior.jls"), Chains) |> DataFrame

# ╔═╡ 04640130-791c-4aac-8af2-450fbcbf7bdb
lsl_trajs = DrWatson.load(datadir("processed", "msl_brick.jld2"), "trajectories")

# ╔═╡ 4450ffd7-679b-49f0-8128-8f2ff2beb906
md"""
## Scenarios

Now we need to create a database of SOWs.
We'll call it 𝐬 (`\bfs`) to mirror the notation of the paper.
"""

# ╔═╡ 5d8adb71-5e3a-41dd-aaad-c76cf3d43c61
n_SOWs = 100_000;

# ╔═╡ fedc0001-84b8-4b2d-856d-34ac7fa99d41
md"""
## Decisions

We'll call our decision space 𝐱, again for consistency.
We have just one decision, which is how high to elevate the house.
We'll discretize this and look at results over all the discrete values for simplicity.
"""

# ╔═╡ 37a5ab0d-db74-4bff-aa12-f79a6404368a
𝐱 = (0:0.5:16)u"ft";

# ╔═╡ 26eb91da-a71b-4eb9-9a8c-a9a336b34446
md"""
## System Model
"""

# ╔═╡ cf4c37c4-0755-4af3-857c-529cce25583d
md"""
### House
We'll start by defining a house.
We'll use the same house for all our calculations.
A house is characterized by its floor plan, in square feet, its dollar value, and its initial elevation relative to the gauge.
Like Zarekarizi et al (2020), we'll consider a house that is below the BFE.
"""

# ╔═╡ fcb1ddc0-e6b8-4847-b994-fce762046f8e
md"""
### Depth-Damage

There are many curves available.
As Zarekarizi et al (2020) describe, this can make a big difference for predicted impacts.
Since our focus here is on sea level rise, we'll stick the HAZUS curve for now -- a proper decision support tool should look for better estimates of depth-damage curves (and, indeed, should consider additional predictors -- see our paper for more discussion).
"""

# ╔═╡ 9a126e57-0316-4caf-ab6a-89e502fb71f1
md"""
### System Function

Following our paper's notation, we need to define a function 𝑓 (`\itf`) that takes in a scenario sᵢ (`s\_i`) and a decision xⱼ (`x\_j`) and spits out an outcome, which is just a vector of length $K$.
We'll use some notation magic with the `Tullio.jl` package.
"""

# ╔═╡ 0fbff6fd-69b8-4aa9-904a-cf6b004428ba
K = 4;

# ╔═╡ e88d50e9-fe5c-46e8-9cd0-45cebdb0d652
function 𝑓(sᵢ, xⱼ)::Vector
	return 1:K
end

# ╔═╡ f180e72b-d2c5-46d7-a2ef-84c012e59cf4
md"""
## Calculate Outcomes
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
	first(𝐬)
end

# ╔═╡ 31b3c3af-32d6-410d-bca2-f200e0326801
y = zeros(length(𝐬), length(𝐱), K);

# ╔═╡ b04b1325-df46-47c1-9d75-0ce91434c600
y[1, 1, :]

# ╔═╡ 58e0aea3-3ea2-4584-a70e-c670a7184395
@tullio y[i, j, :] = 𝑓(𝐬[i], 𝐱[j]);

# ╔═╡ 5f18caa2-9916-4d24-ac6e-0788357d341d
md"""
Parse the file containing the parameters of the Europa depth-damage model
See [https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R](https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R) for details.
"""

# ╔═╡ a54b43fd-bc71-44bd-a051-98401f8a549d
function parse_europa()
    fname = datadir("raw", "fragility_europa.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = dat[!, "depth_feet"]u"ft"
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end;

# ╔═╡ 417aef6e-f03d-4d8b-998c-b86fb09f403c
md"""
Parse the file containing the parameters of the HAZUS depth-damage model
See [https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R](https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R) for details.
"""

# ╔═╡ 457888c0-8877-4a0f-bb4a-2ece1238b4b2
function parse_hazus()
    fname = datadir("raw", "fragility_hazus.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = (dat[!, "depth_feet"] .+ 0.0)u"ft"
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end;

# ╔═╡ fdf4b10a-bbbb-4ddc-9591-3262afd2d452
md"Generic function to parse the function given a key"

# ╔═╡ 13f2a4a4-39a0-4163-ad0a-a009c41d4fb1
function parse_data(key)
    if key == :europa
        depth, damage_frac = parse_europa()
    elseif key == :hazus
        depth, damage_frac = parse_hazus()
    else
        throw("Invalid key")
    end
    return depth, damage_frac
end;

# ╔═╡ f17cc599-2fff-454a-8597-45ed264829db
md"Get the depth-damage interpolation function"

# ╔═╡ 082e36d3-7c6d-4f26-8bef-a710535272ea
function get_depth_damage(key)
    depth, damage_frac = parse_data(key)
    prepend!(depth, minimum(depth) - 0.1u"ft")
    prepend!(damage_frac, 0)
    interp_fn = Interpolations.LinearInterpolation(
        ustrip.(u"ft", depth),
        damage_frac,
        extrapolation_bc = Interpolations.Flat(),
    )
    damage_fn = function (depth::T) where {T<:Unitful.Length}
        return interp_fn(ustrip.(u"ft", depth))
    end
    return damage_fn
end;

# ╔═╡ 7259cd4f-bdcb-4628-906b-b6c05781ac25
md"All information needed to compute costs for a particular house"

# ╔═╡ 83924484-c6b2-46d1-afc3-78f1b1a0e77b
mutable struct HouseStructure{P<:Unitful.Area,T<:Real,L<:Unitful.Length}
    A::P # the house area, in square feet
    V::T # the structural value of the house, in USD2020
    h::L # the height of the house relative to the gauge
end;

# ╔═╡ 52d492ca-882d-460b-975f-03de1cb273ab
house = HouseStructure(1000u"ft^2", 185_000, 5u"ft");

# ╔═╡ ce166d96-c511-44ad-8758-5e92770a316d
md"copy a house structure"

# ╔═╡ f69af79d-02c3-4dd2-9ab5-595e34eff20b
function copy(hs::HouseStructure)
    return HouseStructure(hs.A, hs.V, hs.h)
end;

# ╔═╡ e9fc3a08-048e-41bf-bd72-f9fa81663cb3
md"""
    elevation_cost(house, Δh)
Compute the cost of elevating a particular HouseStructure `house` by `Δh`, where `house` is a `HouseStructure` and `Δh` is a length.
This cost function follows Zarekarizi et al (2020), which is in turn based on the CLARA model. The valid domain of `Δh` is [0, 14ft].
> Zarekarizi, M., Srikrishnan, V., & Keller, K. (2020). Neglecting uncertainties biases house-elevation decisions to manage riverine flood risks. Nature Communications, 11(1), 5361. https://doi.org/10.1038/s41467-020-19188-9
"""

# ╔═╡ 332ffd8c-f240-4544-b9ea-2eca882200d8
begin
	elevation_thresholds = [0 5.0 8.5 12.0 14.0][:] # cost is piecewise linear
	elevation_rates = [80.36 82.5 86.25 103.75 113.75][:] # cost per unit increase per area
	elevation_itp = LinearInterpolation(elevation_thresholds, elevation_rates)
end;

# ╔═╡ 04bc3669-07d4-4000-8f91-971bea160cbf
function elevation_cost(house::HouseStructure, Δh::T) where {T<:Unitful.Length}
    area_ft2 = to_sq_feet(house.A)
    base_cost = (10000 + 300 + 470 + 4300 + 2175 + 3500) # all costs are in dollars
    if Δh < 0.0ft
        throw(DomainError(Δh, "Cannot lower the house"))
    elseif Δh ≈ 0.0ft
        cost = 0.0
    elseif 0.0ft < Δh <= 14.0ft
        rate = elevation_itp(to_feet(Δh))
        cost = base_cost + area_ft2 * rate
    else
        throw(DomainError(Δh, "Cannot elevate >14ft"))
    end
    return cost
end;

# ╔═╡ 01f44422-a92f-42f7-a0c6-86e797f09c3d
md"Depth vs damage as a fraction of house value"

# ╔═╡ 88694dd1-f070-4255-81f0-6379451f01b8
md"Depth vs damage in actual dollars"

# ╔═╡ 87d22dcc-aeaf-4a1c-9b8c-518873b3674d
md"""
The actual callable functions
"""

# ╔═╡ 6bf54483-ec3d-471a-8012-6246e57b8908
begin
	depth_damage_frac_europa = get_depth_damage(:europa)
	depth_damage_frac_hazus = get_depth_damage(:hazus)
end;

# ╔═╡ 23d027d6-5c53-4d0b-997e-eb0db96f2e8f
function depth_damage_frac(
    house::HouseStructure,
    gauge_depth::T,
    model::Symbol,
) where {T<:Unitful.Length}
    flood_depth = gauge_depth - house.h
    if model == :europa
        damage_frac = depth_damage_frac_europa(flood_depth)
    elseif model == :hazus
        damage_frac = depth_damage_frac_hazus(flood_depth)
    else
        throw("Invalid model $model")
    end
    return damage_frac
end;

# ╔═╡ 890ebc2c-89cb-4ce3-92a7-b64f2f656c90
function depth_damage(
    house::HouseStructure,
    gauge_depth::T,
    model::Symbol,
) where {T<:Unitful.Length}
    damage_frac = depth_damage_frac(house, gauge_depth, model)
    return damage_frac * house.V
end;

# ╔═╡ ef2d2c8d-afda-461e-beab-6ce07b8d8650
begin
	function plot_depth_damage()
		depths = (0:0.1:32)u"ft"
		damage = map(d -> depth_damage(house, d, :hazus), depths)
		p = plot(
			xlabel="Gauge Depth [ft]",
			ylabel="Structural Flood Damage [1000 USD]",
			legend=:topleft,
			size=(500, 500),
		)
		plot!(p, ustrip.(u"ft", depths), damage ./ 1_000, label="HAZUS Model")
		return p
	end
	p1 = plot_depth_damage()
	savefig(plotsdir("depth_damage_cost.pdf"))
	p1
end

# ╔═╡ Cell order:
# ╟─49cd8998-378d-11ec-3ceb-811bf9c58697
# ╠═6eaeeaaf-5fe8-44b2-8957-b47918f26b17
# ╟─ab567401-943f-4877-8572-cfa8ffb531d4
# ╠═86eb1175-f48e-456a-a4d4-1ec5d0c0a55e
# ╠═04640130-791c-4aac-8af2-450fbcbf7bdb
# ╟─4450ffd7-679b-49f0-8128-8f2ff2beb906
# ╠═5d8adb71-5e3a-41dd-aaad-c76cf3d43c61
# ╠═89f27639-7b5f-408c-80e7-7cb16cb688c4
# ╟─fedc0001-84b8-4b2d-856d-34ac7fa99d41
# ╠═37a5ab0d-db74-4bff-aa12-f79a6404368a
# ╟─26eb91da-a71b-4eb9-9a8c-a9a336b34446
# ╟─cf4c37c4-0755-4af3-857c-529cce25583d
# ╠═52d492ca-882d-460b-975f-03de1cb273ab
# ╟─fcb1ddc0-e6b8-4847-b994-fce762046f8e
# ╟─ef2d2c8d-afda-461e-beab-6ce07b8d8650
# ╠═9a126e57-0316-4caf-ab6a-89e502fb71f1
# ╠═0fbff6fd-69b8-4aa9-904a-cf6b004428ba
# ╠═e88d50e9-fe5c-46e8-9cd0-45cebdb0d652
# ╠═f180e72b-d2c5-46d7-a2ef-84c012e59cf4
# ╠═31b3c3af-32d6-410d-bca2-f200e0326801
# ╠═58e0aea3-3ea2-4584-a70e-c670a7184395
# ╠═b04b1325-df46-47c1-9d75-0ce91434c600
# ╠═e5456332-766a-48e3-82c9-34224b7589a2
# ╠═16bb254d-3773-4820-b4dd-fe114f7016a4
# ╟─5f18caa2-9916-4d24-ac6e-0788357d341d
# ╠═a54b43fd-bc71-44bd-a051-98401f8a549d
# ╟─417aef6e-f03d-4d8b-998c-b86fb09f403c
# ╠═457888c0-8877-4a0f-bb4a-2ece1238b4b2
# ╠═fdf4b10a-bbbb-4ddc-9591-3262afd2d452
# ╠═13f2a4a4-39a0-4163-ad0a-a009c41d4fb1
# ╟─f17cc599-2fff-454a-8597-45ed264829db
# ╠═082e36d3-7c6d-4f26-8bef-a710535272ea
# ╠═7259cd4f-bdcb-4628-906b-b6c05781ac25
# ╠═83924484-c6b2-46d1-afc3-78f1b1a0e77b
# ╠═ce166d96-c511-44ad-8758-5e92770a316d
# ╠═f69af79d-02c3-4dd2-9ab5-595e34eff20b
# ╠═e9fc3a08-048e-41bf-bd72-f9fa81663cb3
# ╠═332ffd8c-f240-4544-b9ea-2eca882200d8
# ╠═04bc3669-07d4-4000-8f91-971bea160cbf
# ╟─01f44422-a92f-42f7-a0c6-86e797f09c3d
# ╠═23d027d6-5c53-4d0b-997e-eb0db96f2e8f
# ╟─88694dd1-f070-4255-81f0-6379451f01b8
# ╠═890ebc2c-89cb-4ce3-92a7-b64f2f656c90
# ╟─87d22dcc-aeaf-4a1c-9b8c-518873b3674d
# ╠═6bf54483-ec3d-471a-8012-6246e57b8908
