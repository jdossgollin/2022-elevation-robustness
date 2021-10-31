### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 50a61c27-768a-4cf8-8494-e9f99da073c6
begin
	using Pkg
	Pkg.activate("..")
	using Arrow
	using Base: @kwdef
	using ColorSchemes
	using CSV
	using DataFrames
	using Dates
	using Distributions
	using Downloads
	using DrWatson
	using DynamicPPL
	using LinearAlgebra: normalize
	using Optim
	using PlutoUI
	using Plots
	using Plots: mm
	using ProgressBars
	using Random
	using TimeSeries
	using Turing
	using Statistics
	using StatsBase
	using StatsPlots
	using Unitful
	using UnitfulRecipes
	TableOfContents()
end

# ╔═╡ 37b93f78-3750-11ec-2245-4d3c8e4a5b39
md"""
# Notebook 01: Surge Modeling

> Doss-Gollin and Keller, 2021, In Prep. **Do not circulate!**

In this notebook we'll explore the storm surge data and fit an extreme value distribution to it.
"""

# ╔═╡ 57cc380e-c124-488b-af10-a63acab149b0
md"""
## About Pluto

If you've never seen Pluto notebooks before, check out [this video](https://www.youtube.com/watch?v=HiI4jgDyDhY).
Every cell in this notebook (including this Markdown cell!) is produced by code, but sometimes the code is hidden.
To see it, click the eyeball icon to the left.
"""

# ╔═╡ e33670da-89a8-4feb-8b0c-35cead899baf
md"""
## Get data

We're going to work with the time series of annual maximum storm surges, after subtracting historic mean sea level.
We start by reading in the hourly data (see helper functions below for more on how this actually works!)
"""

# ╔═╡ ee7a0e99-3e1b-4413-ac40-0628fa1d5d0b
md"""
We can plot this annual data
"""

# ╔═╡ 140cf28e-e0ce-41f2-a61e-51907805ad40
colors = ColorSchemes.okabe_ito; # colorblind friendly

# ╔═╡ 80915922-3f67-4d2f-853d-6e50c1e8706e
md"""
## Maximum Likelihood Fit

Next, we use the `Extremes` packages to develop a fit to the storm surge data, analogous to the point estimates used in most official guidance.
We'll use the probability-weighted moments method since this is widely used in practice.
"""

# ╔═╡ 806a03f3-13c9-4747-b93c-4a45f796d581
@model function MLEGEV(y)
    μ ~ Turing.Flat()
	ϕ ~ Turing.Flat() # reparameterize
    σ = exp(ϕ)
    ξ ~ Turing.FlatPos(0.0)
    y ~ GeneralizedExtremeValue(μ, σ, ξ)
end;

# ╔═╡ 40b6f48b-efcd-43d0-b6f1-39c40abe3959
md"""
It's always good to look at some diagnostics; these appear qualitatively plausible.
Next we can look at the CDF of this distribution to estimate the 100 year surge.
"""

# ╔═╡ 683e3231-79ea-4f92-93a3-ccac083ec80e
md"""
## Bayesian Fit

In our analysis we want to consider parametric uncertainty in the GEV fit through Monte Carlo sampling.
We will use the Monte Carlo samples to help create our "states of the world (SOWs)".
See [`Turing.jl` docs](https://turing.ml/stable/) for more information on syntax.
"""

# ╔═╡ b1d2dc0d-e456-468e-bda9-c71319ac82bb
md"""
### Model Specification

We use a Generalized Extreme Value (GEV) distribution to model the annual maxima time series.
We fit a stationary model that assumes the data are IID (independent and identically distributed), which is an imperfect assumption.
Future analysis could consider autocorrelation and trends (see something like [Wong, 2018](https://doi.org/10.5194/ascmo-4-53-2018)).

We next need to consider our prior.
Our prior of storm surge will be loosely based on those of Coles & Tawn (1996) and Stephenson (2015), albeit with a few differences.
Specifically:

1. We assume a stationary GEV model for the storm surges $y_i$: $y_i \sim \text{GEV}(\mu, \sigma, \xi)$
1. We assume flat priors on $\mu,\phi := \log \sigma,$ and $\xi$.
1. We impose $\xi>0$ to ensure that the distribution has a lower bound rather than an upper bound. This is reasonable only for our specific circumstances: hourly storm surge plus tide is defined as the residual after subtracting the mean sea level. When we take the maximum of that, it is guaranteed to be $\geq 0$.
1. Following the below references, we impose priors on quantiles of the GEV distribution. That is, rather than try to use domain knowledge to reason about $\mu,\sigma,\xi$ directly, we instead use domain knowledge to reason about the quantiles of the resulting distribution. We apply `LogNormal` priors so that a qiven return level is $\geq 0$ and unbounded.
1. The specific LogNormal parameters are chosen based on weak intuition about storm surges along the US Atlantic and Gulf Coast.

> Coles, S. G., & Tawn, J. A. (1996). A Bayesian analysis of extreme rainfall data. Journal of the Royal Statistical Society: Series C (Applied Statistics), 45(4), 463–478. https://doi.org/10.2307/2986068

> Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC. Retrieved from http://ebookcentral.proquest.com/lib/rice/detail.action?docID=4312572
"""

# ╔═╡ 68e097e1-34f7-4e81-9b60-0c0e7a7631c5
gev_priors = [
	(0.9, Distributions.LogNormal(1.75, 0.2)), # 10 year
	(0.99, Distributions.LogNormal(2, 0.2)), # 100 year
	(0.999, Distributions.LogNormal(2.25, 0.2)), # 1000 year
];

# ╔═╡ a1deb075-4660-4c1c-b084-cc2af869cc21
md"""
>TODO: we can revise these priors
"""

# ╔═╡ b198fc1f-4655-4020-80e6-4d113ed20296
md"""
### Visualize Priors

The prior model described above defines a prior distribution over the 10, 100, and 1000 year floods.
This is an alternative to holding prior beliefs over the parameters of a GEV distribution explicitly.
Let's look at this prior!

Formally, this plot is a *prior predictive check* that shows the distribution of three quantities that depend on the parameters $\theta := (\mu, \sigma, \xi)$
"""

# ╔═╡ 808f2a15-36bb-4855-a044-f86656bfdcfb
begin
	function plot_rl_priors(priors)
		p = plot(
			xlims = (0, 25),
			xlabel = "Storm Surge (ft)",
			ylabel = "Cumulative Probability",
			leftmargin = 5mm,
			legend = :topright,
			yticks = 0:0.2:1,
		)
		for ((prob, dist), c) in zip(priors, colors)
			return_level = Int(round(1 / (1 - prob)))
			plot!(p, x -> pdf(dist, x), label = "$(return_level) Year Flood", linewidth = 2, color=c)
		end
		return p
	end
	p3 = plot_rl_priors(gev_priors)
	savefig(p3, plotsdir("priors_on_surge_quantiles.pdf"))
	plot(p3)
end

# ╔═╡ 600cf274-ff38-42ef-8056-2d3d4de73f2a
md"""
The above quantiles were reached after some iteration and are by no means *correct*, but they seem subjectively defensible.

It's worth noting that we spent quite some time iterating on this model -- we didn't come up with it instantly.
For example, we started with completely uninformative priors and built up incrementally from there.
It probably would have been best for us to retain all these intermediate models (see Gelman et al, 2020) and explain how we reached this model incrementally.

> Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B., Yao, Y., et al. (2020). Bayesian Workflow. ArXiv:2011.01808 [Stat]. Retrieved from http://arxiv.org/abs/2011.01808
"""

# ╔═╡ f76ba3d4-f3f5-4b32-8a84-b4d623c91872
md"""
With that said, here's the model
"""

# ╔═╡ dd688efd-f37a-4dc2-ae39-f5107fd0412e
@model function GEVModel(y)

    # completely flat priors (albeit positive -- see docstring)
    μ ~ Turing.FlatPos(0.0)
    σ ~ Turing.FlatPos(0.0)
    ξ ~ Turing.FlatPos(0.0)

    dist = Distributions.GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for (prob, prior_dist) in gev_priors
        rl = Distributions.quantile(dist, prob)
        Turing.@addlogprob!(logpdf(prior_dist, rl))
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end
end;

# ╔═╡ 72604158-8d54-4cad-a7fa-22fcf3cd47f5
md"""
### Prior Checks
"""

# ╔═╡ 09838a64-d75a-41b8-a518-e7f9fd4c1a50
prior_model = GEVModel(missing);

# ╔═╡ c152a675-097e-4fa6-92f9-714f87c66691
md"We first ask a basic computational question: did the chains mix?"

# ╔═╡ 5956c651-294c-442c-a8ea-ee0635bf6079
md"""
Looks defensible!
It's pretty different from our MLE but the MLE from our data is towards the middle of our prior distribution.
"""

# ╔═╡ 28c54f53-218f-42e4-8e1b-a6c6d7f132ae
md"""
### Synthetic Data Experiment

A desirable attribute of a model of storm surges is being able to accurately fit data that comes from a known 'true' distribution.
We can check that with fake data!
"""

# ╔═╡ 8a33c2e1-d073-438a-a8da-33e3e1f747dc
md"""
This looks consistent with what we know about storm surge in the Eastern United States in general, before looking at data in Norfolk, VA, which is one way to conceptualize what a prior should do.
In theory our model allows negative values of storm surge: although the shape parameter of the GEV is restricted to be $\xi \geq 0$, which ensures that the distribution is left bounded, the lower bound of a GEV distribution is $\mu - \sigma/\xi$, which can be $<0$.
We're not particularly worried about the left tail - the right tail is much more interesting to us.
"""

# ╔═╡ 60e16292-b948-4905-a4fe-605d2bbcf974
md"""
### Fit Posterior

Only now that we have run a few basic checks to ensure that our model seems plausible are we ready to actually look at data!
"""

# ╔═╡ b5812bc5-a1e3-4538-95f7-7dca2758cbef
N_draws = 10_000;

# ╔═╡ 44108282-59c1-477a-84b1-5d0176765476
md"""
Unsurprisingly, our uncertainties have shrunk well beyond the prior.
The biggest uncertainty is in our shape parameter -- this is unsurprising given prior research, which emphasizes that uncertainty in the shape parameter is often large.

We can see this in the predictive distribution for storm surge
"""

# ╔═╡ 441a9767-c82d-4f77-a7b2-7aab95e6fb0f
md"""
### Posterior Test Statistics

Now that we have a posterior distribution, we can compute some test statistics to see how the distribution of these statistics computed over samples from our posterior predictive distribution compare to the values in our observed data.

For a more nuanced discussion, see

> Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(2), 389–402. https://doi.org/10.1111/rssa.12378
"""

# ╔═╡ 1473baa8-98ac-4cef-81c6-1b7044f0e599
md"""
We can see that our estimates (histogram and density) are broadly quite consistent with the observed data (vertical green line). The main discrepancies are:

- We seem to consistently underestimate the (absolute value of) partial autocorrelation, which makes sense: we fit an IID model, which is an approximation of the time dependence of the actual climate. There is some evidence for time dependence in this dataset, but the difference is not huge.
- We have a non-zero (albeit very low!) probability of a storm surge in excess of 30 feet. This is not very realistic. However, since we will be taking the convolution of hazard with a damage function that is S-shaped (even if the house is fully submerged, damage cannot be more than the house value), it doesn't matter *for our purposes here*  whether a flood is 30 or 60 feet, so this is not a fatal flaw
"""

# ╔═╡ 97cb91a4-abde-4752-a13a-343ae69a2d5b
md"""
### Posterior Plots 

Finally we'll create a plot of draws from our posterior
"""

# ╔═╡ 128fccc1-fab6-4a57-9b9c-76ce851ba503
md"""
we can also visualize this as a return period
"""

# ╔═╡ 2ddb440a-812b-4766-b4f1-73e5bb270ead
md"""
Last but not least, we can visualize the spread of the 100 year return level in our posterior
"""

# ╔═╡ 5586c372-498e-4649-a455-b7054ca33404
md"""
## Save Posterior

Now we're all done!
We can save our posterior for easy access in subsequent notebooks
"""

# ╔═╡ acf6d7d6-8ab7-4015-adf9-f4245af9ee54
md"""
## Helper Functions

The functions below help us work with the code we need.
Pluto notebooks do some fancy tricks to let us store this code at the bottom of the file, even though we need to run them before we run some of the other blocks.
"""

# ╔═╡ 12fc0c0d-206d-4b74-b2b5-4d527cc0fadd
md"data structure to hold hourly data"

# ╔═╡ 0ba3cc40-2bdb-45d3-92d0-6a2ae6a7d0c1
struct HourlyGageRecord
    data::TimeSeries.TimeArray
    gage_name::String
    gage_id::String
end

# ╔═╡ 1767abfd-7925-40b0-aac6-ef36c3056965
md"""
Get the URL for the NOAA file containing tides and currents data
"""

# ╔═╡ edca158d-9a45-4572-8d74-b10ab776f34f
function get_noaa_url(year::Int; datum="NAVD", station="8638610", units="metric", time_zone="GMT")
    return "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=hourly_height&application=NOS.COOPS.TAC.WL&begin_date=$(year)0101&end_date=$(year)1231&datum=$(datum)&station=$(station)&time_zone=$(time_zone)&units=$(units)&format=csv"
end;

# ╔═╡ 5970d43e-ff14-4fac-9827-8ddc330d2b60
md"Get the path of local CSV file where to save a single year of data"

# ╔═╡ b0cedf53-fdb7-47e0-8e96-5ce5bccce2e9
function get_csv_fname(year::Int)
    cache_dir = datadir("processed", "gauge")
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    return joinpath(cache_dir, "$(year).csv")
end;

# ╔═╡ 3e55612d-8c26-44cb-b566-ff0891a86107
md"Parse one of the CSV files downloaded from NOAA"

# ╔═╡ 0361bed0-96a9-4c96-964d-4778a88fea77
function parse_noaa_csv_file(fname)
    surge = DataFrame(CSV.File(fname))[!, [1, 2]]
    DataFrames.rename!(surge, [:datetime, :sl])
    surge[!, :datetime] =
        map(x -> DateTime(x, "YYYY-mm-dd HH:MM"), surge[!, :datetime])
    surge[!, :sl] = surge[!, :sl]u"m" # explicitly specify unit to avoid issues
    surge[!, :sl] = uconvert.(Unitful.u"ft", surge[!, :sl])
    dropmissing!(surge)
    return surge
end;

# ╔═╡ 621a55b6-9d0f-4fb5-bc11-cbe1bdf74a13
md"Get a DataFrame of hourly sea level readings"

# ╔═╡ 27242ec1-c969-4416-b55d-50d07a15b655
function get_hourly_obs(year::Int)
    fname = get_csv_fname(year)
    try
        ts = parse_noaa_csv_file(fname)
        return ts
    catch
        url = get_noaa_url(year)
        @info "downloading data for $year"
        Downloads.download(url, fname)
        ts = parse_noaa_csv_file(fname)
        return ts
    end
end;

# ╔═╡ 5945cda0-53a7-4a44-b5a7-9c153e900c75
md"Get all hourly sea level data, spefically for the Sewell's Point, VA (Norfolk) data"

# ╔═╡ 1d0adce9-877f-412f-9b5e-8a9112603c2b
function produce_norfolk_hourly()
    years = 1928:2020
    df = vcat([get_hourly_obs(year) for year in ProgressBars.tqdm(years)]...)
    hourly = TimeSeries.TimeArray(df, timestamp = :datetime)[:sl]
    return HourlyGageRecord(hourly, "Sewell's Point, VA", "8638610")
end;

# ╔═╡ 53444f1f-51f0-4557-871f-c822890659b6
md"convenience function that will save the hourly data"

# ╔═╡ 3db5ff1a-cbe3-4c9f-9ef8-e715691344b8
function get_norfolk_hourly(; overwrite::Bool = false)

    cache_dir = datadir("processed")
    cachename = joinpath(cache_dir, "hourly_surge.jld2")

    try
        @assert !overwrite
        hourly = DrWatson.load(cachename, "hourly")
        return hourly
    catch err
        hourly = produce_norfolk_hourly()
        DrWatson.wsave(cachename, Dict("hourly" => hourly))
        return hourly
    end
end;

# ╔═╡ 4b1871df-02f7-4604-9751-5909d4d8bbe1
hourly = get_norfolk_hourly()

# ╔═╡ 6dc91ac0-94e8-4461-86a9-8a593283bd80
md"Create a structure for annual gauge records"

# ╔═╡ d9f1091d-3209-4b6f-a9fc-0eada1ae881d
@kwdef struct AnnualGageRecord
    year::Vector{<:Int}
    surge::Vector{<:Unitful.Length}
    msl::Vector{<:Unitful.Length}
    gage_name::String
    gage_id::String
end;

# ╔═╡ 0b23d879-f8b7-45be-9a01-b922cc1292d6
md"Convert an hourly to annual record"

# ╔═╡ 6a03d552-d22f-4846-aebf-862230462063
function hourly_to_annual(hourly)
    hourly_df = DataFrame(hourly.data)
    DataFrames.dropmissing!(hourly_df)
    hourly_df[!, :year] = [year(date) for date in hourly_df[!, :timestamp]]
    annual_df = DataFrames.combine(
        DataFrames.groupby(hourly_df, :year),
        :sl => maximum => :surge,
        :sl => mean => :msl,
    )
    annual_df[!, :gage_id] .= hourly.gage_id
    annual_df[!, :gage_name] .= hourly.gage_name
    annual_df[!, :surge] = annual_df[!, :surge] - annual_df[!, :msl]
    return AnnualGageRecord(
        year = annual_df[!, :year],
        surge = annual_df[!, :surge],
        msl = annual_df[!, :msl],
        gage_name = annual_df[!, :gage_name][1],
        gage_id = annual_df[!, :gage_id][1],
    )
end;

# ╔═╡ f9a2c638-2988-4a70-b774-5dc6babc3aaf
annual = hourly_to_annual(hourly)

# ╔═╡ 669536c8-db51-479c-924f-bbad20faae02
annmax_surge_ft = ustrip.(u"ft", annual.surge);

# ╔═╡ f1743aac-c262-4b51-90b9-bcfe7b1533b2
mle_estimate = optimize(MLEGEV(annmax_surge_ft), MLE(), [1, 1, 0.5])

# ╔═╡ 92859431-1be7-4e46-809a-88ac17c968e2
mle = GeneralizedExtremeValue(
	mle_estimate.values[:μ],
	exp(mle_estimate.values[:ϕ]),
	mle_estimate.values[:ξ],
);

# ╔═╡ b4966d0e-6825-47f1-ba49-5cea35a8e342
n_years = length(annmax_surge_ft);

# ╔═╡ 6d808cfd-10df-4107-8b27-2e9a20bbcdb6
md"""
By looking at this plot, we can see that in this one case, we are able to recover the truth quite nicely with $n_years data points, just like in our observed data.
We note that our prior distribution has introduced a very slight bias, in which the estimates of return levels for short return periods (<50 years) are very slightly over-estimated, and the estimates for higher return periods (>50 years) are slightly under-estimated.
This is to be expected, and one interpretation of Bayesian inference is as a regularizing device that trades off a small bias in exchange for reduced variance in estimates.
In any case, the differences here are quite small, and we can take this as evidence that our estimation procedure is capable of recovering the true distribution from $n_years data points, though this is no absolute guarantee that it will do so with the data set of interest.

Another prior predictive check is to sample directly from *prior predictive* distribution, which is just
$$p(y) = p(y | \theta) p(\theta)$$.
"""

# ╔═╡ a91ae172-171f-41f0-86b3-e6e5c8adb3ec
posterior_model = GEVModel(annmax_surge_ft);

# ╔═╡ c85a5ad9-7a8f-46b6-af8a-7e86eaac2c63
begin 
	bfe = quantile(mle, 0.99)
	biggest = maximum(ustrip.(u"ft", annual.surge))
end;

# ╔═╡ 3fe5866f-e500-4865-9919-25ddd0ff7b3d
begin
	function plot_surge_cdf()
		xplot = 1.5:0.05:11
		yplot = cdf.(mle, xplot)
		p = plot(
			xplot,
			yplot,
			xlabel="Storm Surge [ft]",
			ylabel="Cumulative Probability",
			legend=:right,
			label="Best Fit",
			linewidth=2,
			color=colors[1],
			size=(500, 500),
		)
		plot!(
			p,
			[minimum(xplot), bfe, bfe],
			[0.99, 0.99, 0],
			label="100 Year Surge",
			color=colors[2],
			linewidth=1.5,
			marker=".",
			markerfacecolor=colors[2],
			markeredgecolor=false,
		)
		vline!(p, [maximum(annmax_surge_ft)], color=colors[3], label="Observed Max")
		return p
	end
	p2 = plot_surge_cdf()
	savefig(p2, plotsdir("surge_cdf.pdf"))
	plot(p2)
end

# ╔═╡ b50a65fe-b2b3-489d-b3ac-1d8515f19b25
md"We see that our estimate of the 100 year flood, $(bfe) ft, is similar to the biggest flood observed, $(biggest)ft."

# ╔═╡ 6b69c1f3-5210-4eab-8939-b04e2761da47
md"Get the annual data"

# ╔═╡ 811fac1c-f8e0-41e7-9dfa-01fe8f7b1f28
function get_norfolk_annual()
	return hourly_to_annual(get_norfolk_hourly())
end;

# ╔═╡ d694f959-4946-488a-a505-5df4ef43297d
begin
	struct HistoricSurge
		name::String
		year::Int
		is_tc::Bool
	end
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
			size = (750, 500),
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
	end
	p1 = plot_historic_annmax_flood()
	savefig(p1, plotsdir("historic_surge.pdf"))
	plot(p1)
end

# ╔═╡ b44b62e5-5f7e-4649-b7aa-89a1c6d633a2
md"Cache the sampling outputs to avoid constant re-running"

# ╔═╡ 153f976b-9be7-4789-82e4-10e6d6770789
function get_fits(
    model::DynamicPPL.Model,
    model_name::String,
    n_samples::Int;
    n_chains::Int = 1,
    drop_warmup::Bool = true,
    overwrite::Bool = false,
)
	cachename = datadir(
		"processed",
		"surge_models",
		"stationary_$(model_name)_$(n_samples).jls",
	)
    samples_per_chain = Int(ceil(n_samples / n_chains))

    try
        @assert !overwrite
        chains = read(cachename, chains)
        return chains
    catch err
        chains = sample(
            model,
            NUTS(),
            MCMCThreads(),
            samples_per_chain,
            n_chains;
            drop_warmup = drop_warmup,
        )
		mkpath(dirname(cachename))
        write(cachename, chains)
        return chains
    end
end;

# ╔═╡ 549f9646-f012-440f-94dc-484e0399d717
prior_fits = get_fits(
	prior_model,
	"prior_model",
	100_000;
	n_chains = 4,
);

# ╔═╡ 41e3000a-1ad3-4f67-b6a3-aa96e706969d
summarystats(prior_fits)

# ╔═╡ 1f95e7a0-c403-4e62-bcc1-bcd523c1c0a9
begin
	p_prior = plot(prior_fits)
	savefig(p_prior, plotsdir("surge_prior_convergence.pdf"))
	p_prior
end

# ╔═╡ 8af5a42e-6a5a-4e81-9f16-a0a6af64fb04
begin
	function plot_fake_data_experiment()
		fake_dist = GeneralizedExtremeValue(4, 0.5, 0.15)
		rng1 = Random.MersenneTwister(713) # Houston Area Code
		fake_data = rand(rng1, fake_dist, length(annmax_surge_ft))
		fake_posterior = get_fits(GEVModel(fake_data), "fake_data", 10_000)
		rng2 = Random.MersenneTwister(203) # New Haven
		fake_yhat = vcat(
			rand.(
				rng2,
				GeneralizedExtremeValue.(
					fake_posterior[:μ],
					fake_posterior[:σ],
					fake_posterior[:ξ],
				),
				10,
			)...,
		)
		p = plot(
			xlabel = "True Return Level [ft]",
			ylabel = "Estimated Return Level [ft]",
			aspect_ratio = :equal,
			legend = :topleft,
			size = (500, 500),
		)
		for rt in [2, 5, 10, 25, 50, 100, 250, 500]
			q = 1 - 1 / rt
			xplot = quantile(fake_dist, q)
			yplot = quantile(fake_yhat, q)
			scatter!(
				p,
				[xplot],
				[yplot],
				label = "",
				markersize = 3,
			)
			annotate!(
				p,
				xplot + 0.225,
				yplot,
				("$rt year", 8, :left)
			)
		end
		Plots.abline!(p, 1, 0, label = "1:1 Line", color=:gray)
		xx = 1 .- 1 ./ exp.(range(log(2), log(1000); length=100))
		plot!(
			p,
			quantile.(fake_dist, xx),
			[quantile(fake_yhat, xxi) for xxi in xx],
			label="Sampled Posterior",
			color = colors[1],
		)
		annotate!(
				p,
				[5.5, 8.5],
				[8.5, 5.5],
				[
					text("Over-Estimate", 8, :left, rotation=45, colors[2]),
					text("Under-Estimate", 8, :left, rotation=45, colors[2]),
				]
			)
		xlims!(p, 3.5, 10.5)
		return p
		end
	p4 = plot_fake_data_experiment()
	savefig(p4, plotsdir("surge_fake_data_experiment.pdf"))
	plot(p4)
end

# ╔═╡ 6c1e1afe-2e07-406d-a3a4-4bbec74b2ecc
posterior_fits = get_fits(
	posterior_model,
	"posterior_model",
	N_draws;
	n_chains = 4,
);

# ╔═╡ 716569f7-d706-4c6c-954c-8ce316596382
summarystats(posterior_fits)

# ╔═╡ 213f25b0-006d-4cda-a360-5511f67970e6
plot(posterior_fits)

# ╔═╡ d9d8c5b1-97b9-4a95-b8e3-34e2285c4156
begin
	function plot_posterior_draws()
		
		df = DataFrame(posterior_fits)
		x = 2:0.05:9.25
		
		p = plot(
			xlabel = "Storm Surge [ft]",
			ylabel = "Cumulative Probability",
			legend = :topleft,
			size = (500, 500),
		)
		for i in rand(1:nrow(df), 1_000)
			dist = GeneralizedExtremeValue(df[i, :μ], df[i, :σ], df[i, :ξ])
			plot!(
				p,
				x,
				cdf.(dist, x),
				color = :gray,
				linewidth = 0.1,
				alpha = 0.25,
				label = false,
			)
		end
		plot!(
			p,
			x,
			cdf.(mle, x),
			color = :blue,
			linewidth = 1.5,
			label = "MLE",
		)
		return p
	end
	p5 = plot_posterior_draws()
	savefig(p5, plotsdir("surge_posterior_draws.png"))
	p5
end

# ╔═╡ a3a7d738-36ea-4057-b9c0-f3293e5aa539
begin
	function weibull_plot_pos(y)
		N = length(y)
		ys = sort(y; rev=false) # sorted values of y
		nxp = xp = [r / (N + 1) for r in 1:N] # exceedance probability
		xp = 1 .- nxp
		return xp, ys
	end
	function plot_return_period(obs)
		
		rts = range(1.25, 1000, length=500) # return periods
		aeps = 1 .- 1 ./ rts # annual exceedance probability
		
		xticks = [2, 5, 10, 25, 50, 100, 250, 500, 1000]

		# create the plot, specifying the axes, grid, and ticks
		p = plot(
			xlabel="Return Period [years]",
			ylabel="Return Level [ft]",
			xscale = :log,
			legend = :topleft,
			xticks = (xticks, string.(xticks)),
			dpi=250, # for saving
		)
		
		# start with Bayesian posterior draws
		df = DataFrame(posterior_fits)
		for i in 1:1_000
			row = df[rand(1:nrow(df)), :]
			dist = GeneralizedExtremeValue(row[:μ], row[:σ], row[:ξ])
			plot!(
				p,
				rts,
				quantile.(dist, aeps),
				color = :gray,
				linewidth = 0.1,
				alpha = 0.25,
				label = (i == 1) ? "Posterior Draws" : false,
			)
		end
		
		# next add the observed data (with Weibull plotting position!)
		xp, ys = weibull_plot_pos(obs)
		scatter!(
			p,
			1 ./ xp,
			ys,
			label="Observed",
			color=colors[1],
			alpha = 1,
		)
		
		# then on top add the MLE
		plot!(
			p,
			rts,
			quantile.(mle, aeps),
			color = colors[2],
			linewidth = 2,
			label = "MLE",
		)

		# return the plot
		return p
	end
	p_rt = plot_return_period(annmax_surge_ft)
	savefig(p_rt, plotsdir("surge_return_period.png"))
	p_rt
end

# ╔═╡ 61a65bef-9e9b-4b00-8c00-1e735b847237
begin
	function plot_100()
		df = DataFrame(posterior_fits)
		rt_100 = map(eachrow(df)) do row
			gev = GeneralizedExtremeValue(row[:μ], row[:σ], row[:ξ])
			quantile(gev, 0.99)
		end
		p = plot(
			xlabel = "100 Year Return Level [ft]",
			ylabel = "Probability Density",
			size = (500, 300),
		)
		density!(
			p,
			rt_100,
			c=colors[2],
			label="Posterior Draws",
			linewidth=3,
		)
		vline!(p, [quantile(mle, 0.99)], c=colors[1], label="MLE", linewidth=3)
		vline!(
			p,
			[maximum(annmax_surge_ft)],
			c=colors[3],
			label="Historical Record",
			linewidth=3,
		)
		return p
	end
	p_100 = plot_100()
	savefig(p_100, plotsdir("surge_100_year_uncertainty.pdf"))
	p_100
end

# ╔═╡ d0fe9336-6eb6-404c-bcd3-cc5299b1c464
write(
	datadir("processed", "surge_posterior.jls"),
	posterior_fits
)

# ╔═╡ 59213bea-dd9c-45bd-84b0-a1fd731a0452
md"Draws samples from the posterior distribution of fit. This will simulate a vector of `N` samples for each  posterior draw in `fit`."

# ╔═╡ 1e2137bf-c514-4948-886e-c3f9163b25ee
function sample_predictive(fit::Chains, N::Int)
	df = DataFrame(fit)
    yhat = [
        rand(GeneralizedExtremeValue(μ, σ, ξ), N) for
        (μ, σ, ξ) in zip(df[!, :μ], df[!, :σ], df[!, :ξ])
    ]
    return yhat
end;

# ╔═╡ e5ac70c4-56fc-4361-84cb-34098e456053
begin
	function plot_predictive(p, chains::Chains, N::Int, label::String)
		predictive = sample_predictive(chains, N)
		predictive_flat = vcat(predictive...)[:]
		edges = range(
			quantile(predictive_flat, 0.001),
			quantile(predictive_flat, 0.999),
			length=250,
		)
		h = normalize(fit(Histogram, predictive_flat, edges); mode=:pdf)
		r = h.edges[1]
		x = first(r)+step(r)/2:step(r):last(r)
		plot!(
			p,
			x,
			h.weights,
			label=label,
			xlabel="Storm Surge [ft]",
			ylabel="Probability Density",
			left_margin=2.5mm,
		)
		return p
	end
	function plot_predictive(chains::Chains, N::Int, label::String)
		p = plot()
		return plot_predictive(p, chains, N, label)
	end
	plot_predictive(prior_fits, n_years, "prior predictive")
end

# ╔═╡ b3ddd371-7d1a-48a9-b88d-c8ae34f3a5a2
begin
	p_predictive = plot_predictive(prior_fits, n_years, "Prior Predictive")
	plot_predictive(p_predictive, posterior_fits, n_years, "Posterior Predictive")
	savefig(p_predictive, plotsdir("prior_posterior_predictive.pdf"))
	p_predictive
end

# ╔═╡ 0f93af74-6cf8-4e15-a349-3c137c109b3a
begin
	function plot_test_stat(t, yhat; title="", xlabel=:"")
		t_posterior = [t(yi) for yi in yhat]
		observed = t(annmax_surge_ft)
		bins=range(
			quantile(t_posterior, 0.001),
			quantile(t_posterior, 0.999),
			length=100,
		)
		p = histogram(
			t_posterior,
			bins = bins,
			title = title,
			xlabel = xlabel,
			label = "PPD",
			normalize = :pdf,
			yticks = :none,
			linealpha=0,
			fillcolor = colors[1],
		)
		vline!(p, [observed], linewidth = 4, label = "obs", color = colors[2])
		return p
	end
	function MannKendall(x::Vector{<:Real})
		return corkendall([1:length(x);], x)
	end
	test_stats = [
		("Lag 1 PACF", "Autocorrelation", x -> pacf(x, 1:1)[1]),
		("Lag 2 PACF", "Autocorrelation", x -> pacf(x, 2:2)[1]),
		("Lag 3 PACF", "Autocorrelation", x -> pacf(x, 3:3)[1]),
		("Lag 5 PACF", "Autocorrelation", x -> pacf(x, 5:5)[1]),
		("Maximum Flood", "Surge (ft)", maximum),
		("Minimum Flood", "Surge (ft)", minimum),
		("Median Flood", "Surge (ft)", median),
		("Third Biggest Flood", "Surge (ft)", x -> reverse(x[sortperm(x)])[3]),
		("Mann-Kendall Test", "Test Statistic", MannKendall),
	]
	yhat = sample_predictive(posterior_fits, n_years)
	test_plots = plot(
		[
			plot_test_stat(t, yhat, title=title, xlabel=xlabel)
			for (title, xlabel, t) in test_stats
		]...,
		size = (1200, 900),
	)
	savefig(test_plots, plotsdir("surge_test_statistics.pdf"))
	plot(test_plots)
end

# ╔═╡ Cell order:
# ╟─37b93f78-3750-11ec-2245-4d3c8e4a5b39
# ╠═50a61c27-768a-4cf8-8494-e9f99da073c6
# ╟─57cc380e-c124-488b-af10-a63acab149b0
# ╟─e33670da-89a8-4feb-8b0c-35cead899baf
# ╠═4b1871df-02f7-4604-9751-5909d4d8bbe1
# ╠═f9a2c638-2988-4a70-b774-5dc6babc3aaf
# ╟─ee7a0e99-3e1b-4413-ac40-0628fa1d5d0b
# ╠═140cf28e-e0ce-41f2-a61e-51907805ad40
# ╟─d694f959-4946-488a-a505-5df4ef43297d
# ╟─80915922-3f67-4d2f-853d-6e50c1e8706e
# ╠═669536c8-db51-479c-924f-bbad20faae02
# ╠═806a03f3-13c9-4747-b93c-4a45f796d581
# ╠═f1743aac-c262-4b51-90b9-bcfe7b1533b2
# ╠═92859431-1be7-4e46-809a-88ac17c968e2
# ╟─40b6f48b-efcd-43d0-b6f1-39c40abe3959
# ╟─c85a5ad9-7a8f-46b6-af8a-7e86eaac2c63
# ╟─3fe5866f-e500-4865-9919-25ddd0ff7b3d
# ╟─b50a65fe-b2b3-489d-b3ac-1d8515f19b25
# ╟─683e3231-79ea-4f92-93a3-ccac083ec80e
# ╟─b1d2dc0d-e456-468e-bda9-c71319ac82bb
# ╠═68e097e1-34f7-4e81-9b60-0c0e7a7631c5
# ╟─a1deb075-4660-4c1c-b084-cc2af869cc21
# ╟─b198fc1f-4655-4020-80e6-4d113ed20296
# ╟─808f2a15-36bb-4855-a044-f86656bfdcfb
# ╟─600cf274-ff38-42ef-8056-2d3d4de73f2a
# ╟─f76ba3d4-f3f5-4b32-8a84-b4d623c91872
# ╠═dd688efd-f37a-4dc2-ae39-f5107fd0412e
# ╟─72604158-8d54-4cad-a7fa-22fcf3cd47f5
# ╠═09838a64-d75a-41b8-a518-e7f9fd4c1a50
# ╠═549f9646-f012-440f-94dc-484e0399d717
# ╠═41e3000a-1ad3-4f67-b6a3-aa96e706969d
# ╟─c152a675-097e-4fa6-92f9-714f87c66691
# ╟─1f95e7a0-c403-4e62-bcc1-bcd523c1c0a9
# ╟─5956c651-294c-442c-a8ea-ee0635bf6079
# ╟─28c54f53-218f-42e4-8e1b-a6c6d7f132ae
# ╟─8af5a42e-6a5a-4e81-9f16-a0a6af64fb04
# ╟─b4966d0e-6825-47f1-ba49-5cea35a8e342
# ╟─6d808cfd-10df-4107-8b27-2e9a20bbcdb6
# ╟─e5ac70c4-56fc-4361-84cb-34098e456053
# ╟─8a33c2e1-d073-438a-a8da-33e3e1f747dc
# ╟─60e16292-b948-4905-a4fe-605d2bbcf974
# ╠═a91ae172-171f-41f0-86b3-e6e5c8adb3ec
# ╠═b5812bc5-a1e3-4538-95f7-7dca2758cbef
# ╠═6c1e1afe-2e07-406d-a3a4-4bbec74b2ecc
# ╟─716569f7-d706-4c6c-954c-8ce316596382
# ╟─213f25b0-006d-4cda-a360-5511f67970e6
# ╟─44108282-59c1-477a-84b1-5d0176765476
# ╟─b3ddd371-7d1a-48a9-b88d-c8ae34f3a5a2
# ╟─441a9767-c82d-4f77-a7b2-7aab95e6fb0f
# ╟─0f93af74-6cf8-4e15-a349-3c137c109b3a
# ╟─1473baa8-98ac-4cef-81c6-1b7044f0e599
# ╟─97cb91a4-abde-4752-a13a-343ae69a2d5b
# ╟─d9d8c5b1-97b9-4a95-b8e3-34e2285c4156
# ╟─128fccc1-fab6-4a57-9b9c-76ce851ba503
# ╟─a3a7d738-36ea-4057-b9c0-f3293e5aa539
# ╟─2ddb440a-812b-4766-b4f1-73e5bb270ead
# ╟─61a65bef-9e9b-4b00-8c00-1e735b847237
# ╟─5586c372-498e-4649-a455-b7054ca33404
# ╠═d0fe9336-6eb6-404c-bcd3-cc5299b1c464
# ╟─acf6d7d6-8ab7-4015-adf9-f4245af9ee54
# ╟─12fc0c0d-206d-4b74-b2b5-4d527cc0fadd
# ╠═0ba3cc40-2bdb-45d3-92d0-6a2ae6a7d0c1
# ╟─1767abfd-7925-40b0-aac6-ef36c3056965
# ╠═edca158d-9a45-4572-8d74-b10ab776f34f
# ╟─5970d43e-ff14-4fac-9827-8ddc330d2b60
# ╠═b0cedf53-fdb7-47e0-8e96-5ce5bccce2e9
# ╟─3e55612d-8c26-44cb-b566-ff0891a86107
# ╠═0361bed0-96a9-4c96-964d-4778a88fea77
# ╟─621a55b6-9d0f-4fb5-bc11-cbe1bdf74a13
# ╠═27242ec1-c969-4416-b55d-50d07a15b655
# ╟─5945cda0-53a7-4a44-b5a7-9c153e900c75
# ╠═1d0adce9-877f-412f-9b5e-8a9112603c2b
# ╟─53444f1f-51f0-4557-871f-c822890659b6
# ╠═3db5ff1a-cbe3-4c9f-9ef8-e715691344b8
# ╟─6dc91ac0-94e8-4461-86a9-8a593283bd80
# ╠═d9f1091d-3209-4b6f-a9fc-0eada1ae881d
# ╟─0b23d879-f8b7-45be-9a01-b922cc1292d6
# ╠═6a03d552-d22f-4846-aebf-862230462063
# ╟─6b69c1f3-5210-4eab-8939-b04e2761da47
# ╠═811fac1c-f8e0-41e7-9dfa-01fe8f7b1f28
# ╟─b44b62e5-5f7e-4649-b7aa-89a1c6d633a2
# ╠═153f976b-9be7-4789-82e4-10e6d6770789
# ╟─59213bea-dd9c-45bd-84b0-a1fd731a0452
# ╠═1e2137bf-c514-4948-886e-c3f9163b25ee
