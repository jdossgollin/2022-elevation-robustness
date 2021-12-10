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
	using Distributions
	using DrWatson
	using Interpolations
	using MCMCChains: Chains
	using Unitful
	using UnitfulRecipes
	using PlutoUI
	using Plots
	using Plots: mm
	using ProgressBars
	using Tullio
	using StaticArrays
	using StatsPlots
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
𝐱 = (0:0.125:14)u"ft";

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
We'll consider a house that was about 1ft below the BFE in 2000.
"""

# ╔═╡ fcb1ddc0-e6b8-4847-b994-fce762046f8e
md"""
### Depth-Damage

There are many curves available.
As Zarekarizi et al (2020) describe, this can make a big difference for predicted impacts.
Since our focus here is on sea level rise, we'll stick the HAZUS curve for now -- a proper decision support tool should look for better estimates of depth-damage curves (and, indeed, should consider additional predictors -- see our paper for more discussion).
"""

# ╔═╡ ba716c60-0b2e-4abe-b195-afa49dbe63d8
md"""
### Elevation Cost

The cost of construction is modeled as piecewise linear as in Zarekarizi et al (2020).
There is no cost of not elevating the house.
The house cannot be lifted more than 14ft.
"""

# ╔═╡ 031e6664-ff59-414b-83eb-1f8eda4c23f7
md"""
### Expected Damage Emulator

We use a Monte Carlo approximation to estimate annual expected damages given a house's elevation relative to mean sea level.

The expected damages ``d`` as a function of ``h``, the difference between the elevation of the house and mean sea level, is the convolution of the probability density of storm surge ``p(y')`` and the depth-damage function ``d(y')``:

``\mathbb{E}(d | h) = \int p(y') d(y') dy'``

This can be approximated as

``
\frac{1}{N} \sum_{i=1}^N d(y^*_i)
``

where ``y^*_1, \ldots, y^*_N`` are draws from ``p(y')``.
With sufficiently large `N`, this is a good approximation.

This is implemented below (see `Helper Functions`).
The function takes in the distance between the mean sea level and the house and returns the expected fraction of damages.
We run this Monte Carlo for a lot of iterations, so it takes quite a while the first time, but then it's cached (we use a progress bar to track... progress... but it shows up in the terminal you're running Pluto from, not on this page. Go check there!).
"""

# ╔═╡ 9a126e57-0316-4caf-ab6a-89e502fb71f1
md"""
### System Function

Following our paper's notation, we need to define a function 𝑓 (`\itf`) that takes in a scenario s and a decision x and spits out an outcome, which is just a vector of length $K$.
"""

# ╔═╡ 5e9e8516-cf5d-4a62-8617-ba7bca1c6580
K = 4; # number of outcomes to track

# ╔═╡ f180e72b-d2c5-46d7-a2ef-84c012e59cf4
md"""
## Calculate Outcomes

We'll use some notation magic with the `Tullio.jl` package to create an array of outcomes $u$ indexed by [scenario, decision, outcome variable].
"""

# ╔═╡ 0bd6b397-59ad-429d-bc84-c5ab6bc91ba2
md"""
## Scenario visualizations
"""

# ╔═╡ adefc367-4f7e-41b1-8666-b716c346cbd1
md"""
## Save Scenarios
"""

# ╔═╡ e5456332-766a-48e3-82c9-34224b7589a2
md"""
## Helper Functions
"""

# ╔═╡ 3e970722-0b5c-4288-80b7-302398d34313
md"""
### Core

Define what a state of the world (``s \in \mathcal{S}``) is
"""

# ╔═╡ 16bb254d-3773-4820-b4dd-fe114f7016a4
struct SOW{I<:Integer,T<:AbstractFloat,L<:Unitful.Length}
    years::UnitRange{I}
    lsl::Vector{L}
    rcp::T
    dynamics::AbstractString
end

# ╔═╡ 89f27639-7b5f-408c-80e7-7cb16cb688c4
begin
	𝐬 = []
	for _ in 1:n_SOWs
		lsl = rand(lsl_trajs)
		sow = SOW(
			lsl.years,
			lsl.lsl,
			lsl.rcp,
			lsl.dynamics,
		)
		push!(𝐬, sow)
	end
	first(𝐬)
end

# ╔═╡ 86277361-9319-46fd-8fe9-448176668e6a
md"""
### Depth Damage
"""

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
md"Generic function to parse the europa/hazus data given a key"

# ╔═╡ f17cc599-2fff-454a-8597-45ed264829db
md"Get the depth-damage interpolation function"

# ╔═╡ 082e36d3-7c6d-4f26-8bef-a710535272ea
function get_depth_damage(key)

	# parse the raw file
	if key == :europa
        depth, damage_frac = parse_europa()
    elseif key == :hazus
        depth, damage_frac = parse_hazus()
    else
        throw("Invalid key")
    end

	# add a zero to the beginning of the data
    prepend!(depth, minimum(depth) - 0.1u"ft")
    prepend!(damage_frac, 0)

	# interpolate
    interp_fn = Interpolations.LinearInterpolation(
        ustrip.(u"ft", depth),
        damage_frac,
        extrapolation_bc = Interpolations.Flat(),
    )

	# return a *function*
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

# ╔═╡ a0a7288a-a838-4f05-950a-68088258d7db
function get_default_house()
	return HouseStructure(1_000u"ft^2", 185_000.0, 7.5u"ft")
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
		damage = map(d -> depth_damage(get_default_house(), d, :hazus), depths)
		p = plot(
			xlabel="Gauge Depth [ft]",
			ylabel="Structural Flood Damage [1000 USD]",
			legend=:topleft,
			size=(400, 400),
		)
		plot!(
			p,
			ustrip.(u"ft", depths),
			damage ./ 1_000,
			label="HAZUS Model",
			linewidth=2,
		)
		return p
	end
	p1 = plot_depth_damage()
	savefig(plotsdir("depth_damage_cost.pdf"))
	p1
end

# ╔═╡ a52ea94f-3f4a-4a0c-b5ba-56f1600da7e5
md"""

### Elevation Cost

Compute the cost of elevating a particular HouseStructure `house` by `Δh`, where `house` is a `HouseStructure` and `Δh` is a length.
This cost function follows Zarekarizi et al (2020), which is in turn based on the CLARA model. The valid domain of `Δh` is [0, 14ft].
"""

# ╔═╡ bba1a940-a3f2-43df-ba0d-a05c03a829be
begin
	# constants
	elevation_thresholds = [0 5.0 8.5 12.0 14.0][:] # piecewise linear
	elevation_rates = [80.36 82.5 86.25 103.75 113.75][:] # cost /ft / ft^2

	# interpolation
	elevation_itp = LinearInterpolation(elevation_thresholds, elevation_rates)

	# user-facing function
	function elevation_cost(house::HouseStructure, Δh::T) where {T<:Unitful.Length}
	    area_ft2 = ustrip(u"ft^2", house.A)
	    base_cost = (10000 + 300 + 470 + 4300 + 2175 + 3500) # in USD
	    if Δh < 0.0u"ft"
	        throw(DomainError(Δh, "Cannot lower the house"))
	    elseif Δh ≈ 0.0u"ft"
	        cost = 0.0
		elseif 0.0u"ft" < Δh <= 14.0u"ft"
	        rate = elevation_itp(ustrip(u"ft", Δh))
	        cost = base_cost + area_ft2 * rate
	    else
	        throw(DomainError(Δh, "Cannot elevate >14ft"))
	    end
	    return cost
	end
end;

# ╔═╡ 9f8016bd-c104-4d0c-9e5f-61a495f4fdbb
begin
	function plot_elevation()
		Δh = (0:0.5:14)u"ft"
		cost = map(Δhᵢ -> elevation_cost(get_default_house(), Δhᵢ), Δh)
		p = plot(
			xlabel="Height Increase [ft]",
			ylabel="Construction Cost [1000 USD]",
			legend=:topleft,
			size=(500, 500),
		)
		plot!(p, ustrip.(u"ft", Δh), cost ./ 1_000, label=false)
		return p
	end
	p2 = plot_elevation()
	savefig(plotsdir("elevation_cost.pdf"))
	p2
end

# ╔═╡ ea3ec96a-b53e-48f9-82f2-8b615772c110
md"""
### Expected Damage Emulator
"""

# ╔═╡ af851075-eeda-4c9b-a442-5052d171a9ea
md"""We've seen this function before to sample from the predictive distribution"""

# ╔═╡ 0e6125a2-bc5f-4bff-b67f-7f382a4e122b
function sample_predictive(fit::Chains, N::Int)
	df = DataFrame(fit)
    yhat = [
        rand(GeneralizedExtremeValue(μ, σ, ξ), N) for
        (μ, σ, ξ) in zip(df[!, :μ], df[!, :σ], df[!, :ξ])
    ]
    return yhat
end;

# ╔═╡ be8ee31b-fa31-46b8-8665-ca73c9e4e527
md"Get a bunch of synthetic storm surges!"

# ╔═╡ e43561b6-41c8-41f0-8d1d-8f5133e9a2f2
function draw_surges(N::Int)
    posterior = read(datadir("processed", "surge_posterior.jls"), Chains)
    y_hat = sample_predictive(posterior, N)
    return vcat(y_hat...) .* 1.0u"ft"
end;

# ╔═╡ 6d90dcaa-a438-49d3-bf22-97a7d897479d
md"""
Fit an interpolation function that estimates expected annual damage given difference between house elevation and MSL.
This returns a function that takes in the distance between the house floor and the mean sea level, and returns an expected annual damage fraction.
The nice thing is that nothinga about this emulator is specific to a particular house (just to a particular depth-damage function).
"""

# ╔═╡ 93e3cbe4-cf3c-4cd4-b5e6-6d40c6b2b78e
function fit_expected_damage_emulator(key::Symbol=:hazus; N::Int=2_000)
    clearances = collect(-30:0.5:30)u"ft"
    surges = draw_surges(N)
	dmg_fn = get_depth_damage(key)
    xp_dmg =
		[
			mean(dmg_fn.(surges .- cᵢ))
			for cᵢ in ProgressBar(clearances)
		]
    interp_fn = Interpolations.LinearInterpolation(
        ustrip.(u"ft", clearances),
        xp_dmg,
        extrapolation_bc = Interpolations.Flat(),
    )
    damage_fn = function (h::T) where {T<:Unitful.Length}
        return interp_fn(ustrip.(u"ft", h))
    end
    return damage_fn
end;

# ╔═╡ 5f0f47a9-c4b0-49bf-b17e-6c8ba7a7ebd5
md"Get an interpolation function that estimates expected annual damage given difference between house elevation and MSL"

# ╔═╡ 5da48412-b3fa-407d-8c3b-63d118f77735
function get_expected_damage_emulator(key::Symbol; overwrite::Bool = false)

    # where to save the emulator
    cachename = datadir("processed", "expected_depth_damage_$key.jld2")

    try
        @assert !overwrite
        damage_fn = DrWatson.load(cachename, "damage_fn")
		@assert 0 <= damage_fn(2u"ft") <= 1
        return damage_fn
    
	catch err
        damage_fn = fit_expected_damage_emulator(key)
        DrWatson.wsave(cachename, Dict("damage_fn" => damage_fn))
        return damage_fn
    end
end;

# ╔═╡ d388fd82-fb94-4874-9308-0d7be5d9cd17
expected_damage = get_expected_damage_emulator(:hazus);

# ╔═╡ 331ee751-2cf7-4d4c-b188-e0fce7e47709
begin
	function plot_expected_dmg_msl()
		house = get_default_house()
		clearances = (0:0.25:21)u"ft"
		damages = map(h -> expected_damage(h), clearances) * house.V
		yticks = [5, 10, 100, 1000, 1000, 5_000, 25_000]
		p = plot(
			xlabel="House Height Above MSL [ft]",
			ylabel="Expected Annual Damages [USD]",
			legend=:bottomleft,
			yscale = :log,
			yticks = (yticks, string.(yticks)),
			xticks = 0:3:21,
			size = (500, 500),
		)
		plot!(p, ustrip.(u"ft", clearances), damages, label=false)
		return p
	end
	p3 = plot_expected_dmg_msl()
	savefig(p3, plotsdir("expected_damage_msl.pdf"))
	p3
end

# ╔═╡ e88d50e9-fe5c-46e8-9cd0-45cebdb0d652
function 𝑓(
	s,
	x;
	house = get_default_house(),
	discount_rate = 0.01,
	house_lifetime_yrs = 81,
) where A <: Unitful.Area{Float64} where L <: Unitful.Length{Float64}

	# sanity checks
	@assert length(s.years) >= house_lifetime_yrs
	
	# outcome 1: up front cost
	construction_cost = elevation_cost(house, x)
	house.h += x

	# outcome 2: NPV damages each year
	N = length(s.years)
	clearances = house.h .- s.lsl
	damage_fracs = map(x -> expected_damage(x), clearances)
	damage_usd = damage_fracs .* house.V
	Γ = (1 - discount_rate) .^ collect(0:(N-1))
	npv_damages = Γ .* damage_usd
	total_npv = sum(npv_damages[1:house_lifetime_yrs])

	# outcome 3: total cost
	total_cost = construction_cost + total_npv

	# outcome 4: probability of flooding in final year
	final_damage_frac = last(damage_fracs)

	# make sure it's the right length and type
	return SVector{K, Float64}(
		construction_cost,
		total_npv,
		total_cost,
		final_damage_frac,
	)
end;

# ╔═╡ 31b3c3af-32d6-410d-bca2-f200e0326801
begin
	u = zeros(Float64, length(𝐬), length(𝐱), K)
	@tullio u[i, j, :] = 𝑓(𝐬[i], 𝐱[j])
end;

# ╔═╡ dd4b7937-b95a-4a12-b306-ae9bfdcd4c3d
begin
	function plot_scenario_map(x, title)
		j = findfirst(x .== 𝐱)
		k = 3 # total costs
		sidx = rand(1:length(𝐬), 25_000)
		ss = 𝐬[sidx]
		msl_2100_ft = ustrip.(u"ft", [last(s.lsl) for s in ss])
		total_costs = u[sidx, j, k]
		p = plot(
			xlabel="MSL at Norfolk, VA in 2100 [ft]",
			ylabel = "Expected NPV Total Costs [1000 USD]",
			title="$title",
			size=(600, 600),
			left_margin = 5mm,
			title_align = :left,
			bottom_margin = 6mm,
		)

		scatter!(
			p,
			msl_2100_ft,
			total_costs ./ 1_000,
			markerstrokewidth=0,
			markersize=1,
			label=false,
			alpha = 0.75,
		)
		return p
	end
	yticks = collect(0:100:600)
	heights_plot = [0, 2.5, 5, 8, 11, 14]u"ft"
	titles = ["($('a'+(i-1))): Δh = $Δh" for (i, Δh) in enumerate(heights_plot)]
	p_scenario_maps = plot(
		[
			plot_scenario_map(xi, t)
			for (xi, t) in zip(heights_plot, titles)
		]...,
		dpi=150,
		size=(1200, 900),
		link=:y,
		ylims=(0, 600),
		yscale=:log10,
		yticks = (yticks, string.(yticks)),
	)
	savefig(p_scenario_maps, plotsdir("scenario_maps.png"))
	p_scenario_maps
end

# ╔═╡ 3f3a888b-4cc9-4415-9e50-ef941b7fdd95
DrWatson.wsave(datadir("processed", "sows_outcomes.jld2"),
	Dict("u" => u,  "𝐬" => 𝐬, "𝐱" => 𝐱)
);

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
# ╠═a0a7288a-a838-4f05-950a-68088258d7db
# ╟─fcb1ddc0-e6b8-4847-b994-fce762046f8e
# ╠═ef2d2c8d-afda-461e-beab-6ce07b8d8650
# ╟─ba716c60-0b2e-4abe-b195-afa49dbe63d8
# ╟─9f8016bd-c104-4d0c-9e5f-61a495f4fdbb
# ╟─031e6664-ff59-414b-83eb-1f8eda4c23f7
# ╠═d388fd82-fb94-4874-9308-0d7be5d9cd17
# ╟─331ee751-2cf7-4d4c-b188-e0fce7e47709
# ╟─9a126e57-0316-4caf-ab6a-89e502fb71f1
# ╠═5e9e8516-cf5d-4a62-8617-ba7bca1c6580
# ╠═e88d50e9-fe5c-46e8-9cd0-45cebdb0d652
# ╟─f180e72b-d2c5-46d7-a2ef-84c012e59cf4
# ╠═31b3c3af-32d6-410d-bca2-f200e0326801
# ╟─0bd6b397-59ad-429d-bc84-c5ab6bc91ba2
# ╠═dd4b7937-b95a-4a12-b306-ae9bfdcd4c3d
# ╟─adefc367-4f7e-41b1-8666-b716c346cbd1
# ╠═3f3a888b-4cc9-4415-9e50-ef941b7fdd95
# ╟─e5456332-766a-48e3-82c9-34224b7589a2
# ╟─3e970722-0b5c-4288-80b7-302398d34313
# ╠═16bb254d-3773-4820-b4dd-fe114f7016a4
# ╟─86277361-9319-46fd-8fe9-448176668e6a
# ╟─5f18caa2-9916-4d24-ac6e-0788357d341d
# ╠═a54b43fd-bc71-44bd-a051-98401f8a549d
# ╟─417aef6e-f03d-4d8b-998c-b86fb09f403c
# ╠═457888c0-8877-4a0f-bb4a-2ece1238b4b2
# ╟─fdf4b10a-bbbb-4ddc-9591-3262afd2d452
# ╟─f17cc599-2fff-454a-8597-45ed264829db
# ╠═082e36d3-7c6d-4f26-8bef-a710535272ea
# ╟─7259cd4f-bdcb-4628-906b-b6c05781ac25
# ╠═83924484-c6b2-46d1-afc3-78f1b1a0e77b
# ╟─01f44422-a92f-42f7-a0c6-86e797f09c3d
# ╠═23d027d6-5c53-4d0b-997e-eb0db96f2e8f
# ╟─88694dd1-f070-4255-81f0-6379451f01b8
# ╠═890ebc2c-89cb-4ce3-92a7-b64f2f656c90
# ╟─87d22dcc-aeaf-4a1c-9b8c-518873b3674d
# ╠═6bf54483-ec3d-471a-8012-6246e57b8908
# ╟─a52ea94f-3f4a-4a0c-b5ba-56f1600da7e5
# ╠═bba1a940-a3f2-43df-ba0d-a05c03a829be
# ╟─ea3ec96a-b53e-48f9-82f2-8b615772c110
# ╟─af851075-eeda-4c9b-a442-5052d171a9ea
# ╠═0e6125a2-bc5f-4bff-b67f-7f382a4e122b
# ╟─be8ee31b-fa31-46b8-8665-ca73c9e4e527
# ╠═e43561b6-41c8-41f0-8d1d-8f5133e9a2f2
# ╟─6d90dcaa-a438-49d3-bf22-97a7d897479d
# ╠═93e3cbe4-cf3c-4cd4-b5e6-6d40c6b2b78e
# ╟─5f0f47a9-c4b0-49bf-b17e-6c8ba7a7ebd5
# ╠═5da48412-b3fa-407d-8c3b-63d118f77735
