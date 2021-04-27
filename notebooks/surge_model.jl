### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 42f6d411-aa95-4dd2-9b7c-b8c94fdd60c5
using DrWatson

# ╔═╡ 7eb464e9-4e88-4980-8cd1-646fe5403a65
begin
	using PlutoUI
	PlutoUI.TableOfContents(title="Outline")
	
	using ColorSchemes
	using CSV
	using DataFrames
	using Plots
	using StatsPlots
	using Turing
	using Unitful
	
	using StatsBase: pacf, corkendall
	using Printf: @sprintf
	
	using NorfolkFloods
	using NorfolkBRICK
end

# ╔═╡ 89ee305e-a1f1-11eb-2b0b-bf4e4913ac08
md"
# Storm Surge Modeling

In this notebook, we will develop and check our model of storm surge for Norfolk, VA.
"

# ╔═╡ 0f50e342-3d55-49af-b611-a43a03ac9fcf
md"## Front Matter

Let's do our package imports and create a table of contents
"

# ╔═╡ cb71e853-b247-4212-804f-0973f03ae2ab
@quickactivate "2021-elevation-robustness"

# ╔═╡ 3bbdc4be-4ea9-427d-ab03-f41b7060f3ab
md"Set a consistent color scheme"

# ╔═╡ 87661bf8-efd7-432b-a89f-00d0564194aa
colors = ColorSchemes.tab10; # define color scheme

# ╔═╡ 90879815-a293-477d-b52b-bc326ad8a941
md"choose how many samples to create and how many chains to run"

# ╔═╡ 71950ea2-aae3-40a6-b46c-b3f61ca73aaf
n_chains, n_samples = 4, 100_000;

# ╔═╡ 42f2eae8-85eb-44a1-8e97-8cddc0722561
md"## Prior Model

We start by exampining our prior.
This prior is provided in the `NorfolkFloods` module included as part of this repository.
However, the full model is specified below here.

Our model of storm surge will be loosely based on those of Coles & Tawn (1996) and Stephenson (2015), albeit with a few differences.
Specifically:

1. We assume a stationary GEV model for the storm surges $y_i$: $y_i \sim \text{GEV}(\mu, \sigma, \xi)$
1. We assume flat priors on $\mu,\sigma,\xi$, but impose $\sigma>0$ by definition, $\xi>0$ to ensure that the distribution has a lower bound, and $\mu > 0$ to restrict the lower bound. This is reasonable only for our specific circumstances: hourly storm surge plus tide is defined as the residual after subtracting the mean sea level. When we take the maximum of that, it is guaranteed to be $\geq 0$.
1. Following the below references, we impose priors on quantiles of the GEV distribution. That is, rather than try to use domain knowledge to reason about $\mu,\sigma,\xi$ directly, we instead use domain knowledge to reason about the quantiles of the resulting distribution. We apply LogNormal priors so that a qiven return level is $\geq 0$ and can potentially be very large. The specific LogNormal parameters are chosen based on weak intuition about storm surges along the US Atlantic and Gulf Coasta.

- Coles, S. G., & Tawn, J. A. (1996). A Bayesian analysis of extreme rainfall data. Journal of the Royal Statistical Society: Series C (Applied Statistics), 45(4), 463–478. https://doi.org/10.2307/2986068
- Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC. Retrieved from http://ebookcentral.proquest.com/lib/rice/detail.action?docID=4312572
"

# ╔═╡ d821aa4f-fe48-481f-b089-6b424ca76065
GEV_priors = get_GEV_priors();

# ╔═╡ 7423d0c2-cde6-4c6a-820e-26cbb78201ea
md"### Visualize Priors

let's visualize these priors"

# ╔═╡ 79eb6578-2c26-4245-8d19-f3e2b742120e
function plot_rl_priors(priors)
    p = plot(
        xlims = (0, 30),
        xlabel = "Storm Surge (ft)",
        ylabel = "Probability Density",
        title = "Assumed Prior Knowledge",
    )
    for (pr, d) in priors
		rl = Int(round(1 / (1 - pr)))
        plot!(p, d, label = "$(rl) Year Flood", linewidth=2)
    end
    return p
end;

# ╔═╡ 1a666d17-5d3c-4393-bb4c-f8fd92c02461
rl_plot = plot_rl_priors(GEV_priors)

# ╔═╡ 840c1399-c74d-44af-93ce-c5d9286b3155
savefig(rl_plot, plotsdir("priors_on_surge_quantiles.pdf"));

# ╔═╡ aba4a98b-5065-4735-9b41-e3f15ed80a89
md"
The above quantiles were reached after some iteration and are by no means *correct*, but they seem *reasonable* (see Gelman and Shalizi, 2013, for a philosophical discussion of priors).

- Gelman, A., & Shalizi, C. R. (2013). Philosophy and the practice of Bayesian statistics. British Journal of Mathematical and Statistical Psychology, 66(1), 8–38. https://doi.org/10.1111/j.2044-8317.2011.02037.x

### Build Model

Now let's build our prior model and draw some samples from it
"

# ╔═╡ 718a638e-b1a4-4f99-affe-11523181e204
prior_model = GEVModel([missing]);

# ╔═╡ c024becc-77bb-4537-9814-49a3907bc6cb
md"### Sample Model

Let's draw some samples from the prior"

# ╔═╡ 65130fcb-4285-4dbe-b1b1-c3edfc6ca33c
prior = get_fits(prior_model, "prior_model", n_samples; n_chains=n_chains);

# ╔═╡ 665d984c-88fd-4f02-9af7-68492832cfc2
prior

# ╔═╡ 6cf4c506-c934-4fd2-8735-7567548abdf7
md"We can visualize our samples and compute some summary stats as a crude way to check convergence"

# ╔═╡ 27b42c59-530f-42a8-a043-eb99c8b8309e
summarystats(prior)

# ╔═╡ 03ff3279-892a-48f0-8a16-ea0ba34bea69
plot(prior)

# ╔═╡ fe5f4b4b-d317-452c-b3d0-fc8592258183
md"### Fake Data Experiment

A desirable attribute of a model of storm surges, is being able to accurately fit data that comes from a known 'true' distribution.
We can check that with fake data!
"

# ╔═╡ 37508701-6009-4975-b538-08f7a497fc9c
begin
	fake_dist = GeneralizedExtremeValue(4, 0.5, 0.15);
	fake_data = DataFrame(CSV.File("../data/processed/fake_data.csv"))[!, :fake_data];
	fake_posterior = get_fits(GEVModel(fake_data), "fake_data", n_samples);
	fake_yhat = vcat(rand.(GeneralizedExtremeValue.(
				fake_posterior[:μ],
				fake_posterior[:σ],
				fake_posterior[:ξ]), 10)...);
	desired = []
	sampled = []
	for q in 1 .- 1 ./ [2, 5, 10, 25, 50, 100, 250, 500, 1000]
		append!(desired, quantile(fake_dist, q))
		append!(sampled, quantile(fake_yhat, q))
	end
end

# ╔═╡ 87a8ac97-e394-414b-ab86-d0db3e7c84ac
begin
	scatter(sampled, desired, xlabel="Estimated", ylabel="Truth", aspect_ratio=:equal, label="Fake Data Experiment", legend=:topleft)
	Plots.abline!(1, 0, label="1:1 Line")
end

# ╔═╡ 4ef33cff-004b-4b94-bb6e-0b049525f1b6
md"Although we see that we are causing a very slight downward bias in the highest quantiles, the effect looks small so we can conclude that our model *can* recover true hazard
"

# ╔═╡ 9667836a-fd01-455d-9032-1e8cc537b47b
md"### Prior Predictive

Now we will sample from the prior predictive distribution
"

# ╔═╡ babba5ae-ea3d-4912-a363-e4ecbc45c1a7
md"
This looks consistent with what we know about storm surge in the Eastern United States in general, before looking at data in Norfolk, VA.
This distribution technically allows values that are negative (the minimum of the GEV can be negative even if the distribution has a lower bound), but this overall looks plausible.

We will add the posterior predictive distribution to this plot later."

# ╔═╡ 8fa92c19-bbed-43b0-9ef5-245410f28439
md"
## Posterior

Only now that we have run a few basic checks to ensure that our model seems plausible are we ready to actually look at data!

### Data

Let's take a look at our data
"

# ╔═╡ 40c36d4e-be34-4d51-b7f6-1bdb2b42b02d
obs = get_norfolk_annual(); # from custom built NorfolkFloods package

# ╔═╡ 17edbd47-1af7-49a7-a0d9-2ba4144f3f8d
y = ustrip.(u"ft", obs.surge);

# ╔═╡ f462f258-d730-484c-873f-0de883d8edbf
yhat_prior, yhat_prior_flat = sample_predictive(prior, length(y));

# ╔═╡ 28d9b25f-1681-472e-9e2f-d49fc2170d8e
predictive_plot = histogram(
	yhat_prior_flat[findall(-2.5 .< yhat_prior_flat .< 15)],
    label = "Prior",
    normalize = :pdf,
    linecolor = false,
	fill=true,
    fillalpha = 0.5,
    fillcolor = colors[1],
    legend = :bottomright,
	orientation=:horizontal,
	ylims = (-2.5, 15),
	xlims = (0, 0.8),
)

# ╔═╡ 249793f3-137b-4406-9c9b-9fbf24bea2bb
md"let's take a quick look at our data"

# ╔═╡ f9f10ec3-3ccf-4954-9f55-ed40649b88f5
data_plot = plot(
    obs.year,
    y,
    label = "",
    ylabel = "Storm Surge (ft)",
    marker = ".",
)

# ╔═╡ 7353d204-622c-4bf9-8319-9fc0ee2ba470
md"### Fit Model

Let's draw some samples from our posterior"

# ╔═╡ 37977dde-5181-4f95-97fb-f6b80f2da858
posterior = get_fits(GEVModel(y), "posterior", n_samples);

# ╔═╡ 04d8d18b-7e00-436d-82b4-8082539c2aac
md"### Sampling Statistics

Let's do the usual visualizatiuon of our sampler"

# ╔═╡ 305ef8c9-519a-4902-b2b8-0b994d115f6d
summarystats(posterior)

# ╔═╡ 8ce3d09d-df1a-4a29-9d22-67f568b7802b
plot(posterior)

# ╔═╡ 939b5fa3-ef2b-48a2-8ea2-b7ab978b7a14
md"### Posterior Predictive

Now let's look at our posterior predictive distribution
"

# ╔═╡ 0f8d90aa-19b4-47c9-9975-fdb0a63b0cd8
yhat, yhat_flat = sample_predictive(posterior, length(y));

# ╔═╡ 49acc5db-fd4f-4104-806d-8e2dc4b54c36
begin
	histogram!(
		predictive_plot,
		yhat_flat[findall(-2.5 .< yhat_flat .< 20)],
		label = "Posterior",
		normalize = :pdf,
		linecolor = false,
		fill=true,
		fillalpha = 0.5,
		fillcolor = colors[2],
		orientation=:horizontal,
	)
	hline!(predictive_plot, [quantile(yhat_prior_flat, 0.99)], color = colors[1], label = "Prior Q99")
	hline!(predictive_plot, [quantile(yhat_flat, 0.99)], color = colors[2], label = "Posterior Q99")
end

# ╔═╡ a7e84094-bb04-4f68-9c09-52ae3c66c1af
md"### Visualization

let's put these plots together and save them nicely"

# ╔═╡ 49bef039-9439-418f-82b4-7a591580eb82
prior_predictive = plot(
	data_plot,
	predictive_plot,
	link = :y,
	layout = (1, 2),
	size=(600, 400),
	ylims = (-1.25, 10),
)

# ╔═╡ 812dcc5e-22f4-439a-b199-11015c4f142e
md"We can see a notable sharpening from a relatively wide prior to a more narrow posterior."

# ╔═╡ ffab9982-6516-4ca0-8139-2c82c08d398b
savefig(prior_predictive, plotsdir("surge_prior_posterior_dists.pdf"));

# ╔═╡ cff18611-e55d-494a-9c95-0ae07e6f986e
md"## Posterior Predictive Checks

Now that we have a posterior distribution, we can compute some test statistics to see how the distribution of these statistics computed over samples from our posterior predictive distribution compare to the values in our observed data"

# ╔═╡ 9d9bebff-e7a4-4398-a17a-72e18c6ff69a
function plot_test_stat(title, xlabel, t)
    t_posterior = [t(yi) for yi in yhat]
    observed = t(y)
    p = histogram(
        t_posterior,
        title = title,
        xlabel = xlabel,
        label = "PPD",
        normalize = :pdf,
        yticks = :none,
        color = colors[1],
    )
    density!(p, t_posterior, linewidth = 3, label = "", color = colors[2])
    vline!(p, [observed], linewidth = 3, label = "obs", color = colors[3])
    return p
end;

# ╔═╡ 53b7054b-309f-4481-95f0-a20dedf905ab
plot_test_stat(title, t) = plot_test_stat(title, "", t)

# ╔═╡ 4ea2d457-1220-4ffc-9288-7902a1d57105
"""The Mann-Kendall trend test"""
function MannKendall(x::Vector{<:Real})
    return corkendall([1:length(x);], x)
end;

# ╔═╡ 8308dc15-b564-4635-bbda-525804d69815
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
];

# ╔═╡ 690e0fe3-e29b-44f2-83dc-1230b4e3ae81
test_plots = plot(
    [plot_test_stat(title, xlabel, t) for (title, xlabel, t) in test_stats]...,
    size = (1200, 900),
)

# ╔═╡ 3c1ff117-6751-4ff7-bb42-b92519ae94ba
savefig(test_plots,plotsdir("surge_test_statistics.pdf"))

# ╔═╡ 2c9d29f5-d155-4a2d-9d0a-3d92f1775323
md"We can see that our estimates (histogram and density) are broadly quite consistent with the observed data (vertical green line). The main discrepancies are:

- We seem to consistently underestimate the (absolute value of) partial autocorrelation, which makes sense: we fit an IID model, which is an approximation of the time dependence of the actual climate. There is some evidence for time dependence in this dataset, but the difference is not huge.
- We have a non-zero (albeit very low!) probability of a storm surge in excess of 30 feet. This is not very realistic. However, since we will be taking the convolution of hazard with a damage function that is S-shaped (even if the house is fully submerged, damage cannot be more than the house value), it doesn't matter *for our purposes here*  whether a flood is 30 or 60 feet, so this is not a fatal flaw
"

# ╔═╡ 03ab7d5e-96eb-4be3-8794-51f200eb3e63
md"## Return Periods

Many hydrologists are used to seeing our PDFs plotted as return periods.
That's easy enough to do!

We will plot the observations using the simple empirical estimator

$$P = \frac{m}{N+1}$$

and we will plot the posterior estimate from our model
"

# ╔═╡ dc46d7cb-3f21-4bdf-bdc8-e22c13a4fed8
function empirical_return_period(y::Vector{<:Real})
	N = length(y)
	m = [sum(y .> yi) for yi in y]
	rt = 1 ./ ( m ./ (N+1))
	return rt
end;

# ╔═╡ 504ebfe6-b3b6-46a8-a9f2-2b9563370735
function plot_return_period(obs::AnnualGageRecord, posterior::Chains; q = [0.1, 0.9])

    yhat, yhat_flat = sample_predictive(posterior, length(obs.surge))
    return_periods = empirical_return_period(ustrip.(u"ft", obs.surge))
	
    rt_plot = 10 .^ (range(0, log10(500); length = 101)[2:end])
    quantile_plot = 1 .- 1 ./ rt_plot
    rl_plot = [quantile(yhat_flat, pr) for pr in quantile_plot]

    best_fits = hcat(
        [
            quantile(GeneralizedExtremeValue(μ, σ, ξ), quantile_plot) for
            (μ, σ, ξ) in zip(posterior[:μ], posterior[:σ], posterior[:ξ])
        ]...,
    )
    lower = [quantile(best_fits[i, :], minimum(q)) for i = 1:size(best_fits)[1]]
    upper = [quantile(best_fits[i, :], maximum(q)) for i = 1:size(best_fits)[1]]
    ci = Int(round(100 * (maximum(q) - minimum(q))))

    xticks = [1, 5, 10, 25, 50, 100, 250, 500]
    ticklabels = [@sprintf("%d", x) for x in xticks]
    p = plot(
        rt_plot,
        rl_plot,
        formatter = identity,
        xticks = (xticks, ticklabels),
        xaxis = :log,
        legend = :topleft,
        label = "Posterior Fit",
        xlabel = "Return Period (years)",
        ylabel = "Estimated Return Level (ft)",
        title = "Storm Surge at $(obs.gage_name)",
    )
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
    scatter!(p, return_periods, y, label = "Observed")
    return p
end;

# ╔═╡ 464d6b50-f0cc-41fb-af02-7010344d6ea4
rt_plot = plot_return_period(obs, posterior; q = [0.1, 0.9])

# ╔═╡ 4480d8c1-ee98-4148-a506-bd04c45b17c4
savefig(rt_plot, plotsdir("return_level.pdf"));

# ╔═╡ 5a1e2a87-8747-4118-97fc-5ed2d7512f86
md"## Save SOWs

Now that we are relatively happy with our model, we can save our set of synthetic future storm surges to use in our decision model.
"

# ╔═╡ 34645ff8-7e7d-4530-b94c-01c4964d192f
future_surges, _ = sample_predictive(posterior, 100); # 100 years

# ╔═╡ 1be519e0-9a6f-48e7-8e63-41a4730fb96c
future_surges_array = transpose(hcat(future_surges...))

# ╔═╡ 8a367bb7-175d-4702-a0df-be18ab0be84d
md"We can see that the variable we have saved is a 100_000 by 100 matrix. The 100_000 rows are 100_000 samples from the posterior distribution. The 100 columns are 100 future years."

# ╔═╡ 577c21bb-41d1-4d82-a0bb-47df3c2a1502
DrWatson.wsave(datadir("processed", "surge_projections.jld2"), Dict("future_surge" => future_surges_array))

# ╔═╡ Cell order:
# ╟─89ee305e-a1f1-11eb-2b0b-bf4e4913ac08
# ╟─0f50e342-3d55-49af-b611-a43a03ac9fcf
# ╠═42f6d411-aa95-4dd2-9b7c-b8c94fdd60c5
# ╠═cb71e853-b247-4212-804f-0973f03ae2ab
# ╠═7eb464e9-4e88-4980-8cd1-646fe5403a65
# ╟─3bbdc4be-4ea9-427d-ab03-f41b7060f3ab
# ╠═87661bf8-efd7-432b-a89f-00d0564194aa
# ╟─90879815-a293-477d-b52b-bc326ad8a941
# ╠═71950ea2-aae3-40a6-b46c-b3f61ca73aaf
# ╟─42f2eae8-85eb-44a1-8e97-8cddc0722561
# ╠═d821aa4f-fe48-481f-b089-6b424ca76065
# ╟─7423d0c2-cde6-4c6a-820e-26cbb78201ea
# ╠═79eb6578-2c26-4245-8d19-f3e2b742120e
# ╠═1a666d17-5d3c-4393-bb4c-f8fd92c02461
# ╠═840c1399-c74d-44af-93ce-c5d9286b3155
# ╟─aba4a98b-5065-4735-9b41-e3f15ed80a89
# ╠═718a638e-b1a4-4f99-affe-11523181e204
# ╟─c024becc-77bb-4537-9814-49a3907bc6cb
# ╠═65130fcb-4285-4dbe-b1b1-c3edfc6ca33c
# ╠═665d984c-88fd-4f02-9af7-68492832cfc2
# ╟─6cf4c506-c934-4fd2-8735-7567548abdf7
# ╠═27b42c59-530f-42a8-a043-eb99c8b8309e
# ╠═03ff3279-892a-48f0-8a16-ea0ba34bea69
# ╟─fe5f4b4b-d317-452c-b3d0-fc8592258183
# ╠═37508701-6009-4975-b538-08f7a497fc9c
# ╠═87a8ac97-e394-414b-ab86-d0db3e7c84ac
# ╟─4ef33cff-004b-4b94-bb6e-0b049525f1b6
# ╟─9667836a-fd01-455d-9032-1e8cc537b47b
# ╠═f462f258-d730-484c-873f-0de883d8edbf
# ╠═28d9b25f-1681-472e-9e2f-d49fc2170d8e
# ╟─babba5ae-ea3d-4912-a363-e4ecbc45c1a7
# ╟─8fa92c19-bbed-43b0-9ef5-245410f28439
# ╠═40c36d4e-be34-4d51-b7f6-1bdb2b42b02d
# ╠═17edbd47-1af7-49a7-a0d9-2ba4144f3f8d
# ╟─249793f3-137b-4406-9c9b-9fbf24bea2bb
# ╠═f9f10ec3-3ccf-4954-9f55-ed40649b88f5
# ╟─7353d204-622c-4bf9-8319-9fc0ee2ba470
# ╠═37977dde-5181-4f95-97fb-f6b80f2da858
# ╟─04d8d18b-7e00-436d-82b4-8082539c2aac
# ╠═305ef8c9-519a-4902-b2b8-0b994d115f6d
# ╠═8ce3d09d-df1a-4a29-9d22-67f568b7802b
# ╟─939b5fa3-ef2b-48a2-8ea2-b7ab978b7a14
# ╠═0f8d90aa-19b4-47c9-9975-fdb0a63b0cd8
# ╠═49acc5db-fd4f-4104-806d-8e2dc4b54c36
# ╟─a7e84094-bb04-4f68-9c09-52ae3c66c1af
# ╠═49bef039-9439-418f-82b4-7a591580eb82
# ╟─812dcc5e-22f4-439a-b199-11015c4f142e
# ╠═ffab9982-6516-4ca0-8139-2c82c08d398b
# ╟─cff18611-e55d-494a-9c95-0ae07e6f986e
# ╠═9d9bebff-e7a4-4398-a17a-72e18c6ff69a
# ╠═53b7054b-309f-4481-95f0-a20dedf905ab
# ╠═4ea2d457-1220-4ffc-9288-7902a1d57105
# ╠═8308dc15-b564-4635-bbda-525804d69815
# ╠═690e0fe3-e29b-44f2-83dc-1230b4e3ae81
# ╠═3c1ff117-6751-4ff7-bb42-b92519ae94ba
# ╟─2c9d29f5-d155-4a2d-9d0a-3d92f1775323
# ╟─03ab7d5e-96eb-4be3-8794-51f200eb3e63
# ╠═dc46d7cb-3f21-4bdf-bdc8-e22c13a4fed8
# ╠═504ebfe6-b3b6-46a8-a9f2-2b9563370735
# ╠═464d6b50-f0cc-41fb-af02-7010344d6ea4
# ╠═4480d8c1-ee98-4148-a506-bd04c45b17c4
# ╟─5a1e2a87-8747-4118-97fc-5ed2d7512f86
# ╠═34645ff8-7e7d-4530-b94c-01c4964d192f
# ╠═1be519e0-9a6f-48e7-8e63-41a4730fb96c
# ╟─8a367bb7-175d-4702-a0df-be18ab0be84d
# ╠═577c21bb-41d1-4d82-a0bb-47df3c2a1502
