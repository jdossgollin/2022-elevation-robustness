### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ f0b9ecec-e017-11eb-3f88-c72476818de3
begin
    using PlutoUI
    PlutoUI.TableOfContents(title = "Outline")
end

# ╔═╡ 38e43be9-d161-42ef-bd06-5649f52fce97
begin
    using ColorSchemes
    using CSV
    using DataFrames
    using DecisionProblem
	using DrWatson
    using Plots
    using StatsBase
    using StatsPlots
    using Turing
    using Unitful

    using StatsBase: pacf, corkendall
    using Printf: @sprintf
	using Plots: mm

    using NorfolkFloods
    using NorfolkBRICK
end

# ╔═╡ fe66d55d-0565-4740-b0a4-ce82d69cc0e9
md"load some required packages and set up the table of contents"

# ╔═╡ e502336a-d13e-468d-8ff1-b0404e20883e
md"
# Storm Surge Modeling

In this notebook, we will develop and check our model of storm surge for Norfolk, VA.
"

# ╔═╡ b779c4e8-4c80-4958-be3c-28f589555cc3
md"""
## Setup

Let's start with some setup.
First, let's pick a consistent color scheme.
We'll choose one that is colorblind-friendly according to [https://juliagraphics.github.io/ColorSchemes.jl/dev/catalogue/](https://juliagraphics.github.io/ColorSchemes.jl/dev/catalogue/).
"""

# ╔═╡ 760742a2-da3e-4ce4-a565-4f54c36ec005
colors = ColorSchemes.okabe_ito;

# ╔═╡ dd48b4f1-37b5-4698-8b15-bf709ce32dc3
md"Next, we load in raw data and choose how many samples to create and how many chains to run"

# ╔═╡ 56067de6-2d51-48ec-afc3-c42b94642508
obs = NorfolkFloods.get_norfolk_annual(); # observed MSL and surge

# ╔═╡ 69b193fa-0d8a-4c4e-81ad-17a776ba7e68
brick = NorfolkBRICK.get_norfolk_brick(syear = obs.year[1]; eyear = obs.year[1] + 100); # projections

# ╔═╡ 0ffeb315-7b5f-4cf4-986e-daae99eecece
n_chains = 4;

# ╔═╡ 190eb3a5-91ec-4a40-a3f2-f816daee1e8f
n_samples = length(brick);

# ╔═╡ c08683e3-718d-400f-9cfc-564d095a3cac
n_years = length(obs.year);

# ╔═╡ 6d211cc4-1c17-4788-be8d-972bfe7d3335
md"## Prior Model

We start by exampining our prior.
This prior is provided in the `NorfolkFloods` module included as part of this repository.
However, the full model is described here as well.

Our model of storm surge will be loosely based on those of Coles & Tawn (1996) and Stephenson (2015), albeit with a few differences.
Specifically:

1. We assume a stationary GEV model for the storm surges $y_i$: $y_i \sim \text{GEV}(\mu, \sigma, \xi)$
1. We assume flat priors on $\mu,\sigma,\xi$, but impose $\sigma>0$ by definition, $\xi>0$ to ensure that the distribution has a lower bound, and $\mu > 0$ to restrict the lower bound. This is reasonable only for our specific circumstances: hourly storm surge plus tide is defined as the residual after subtracting the mean sea level. When we take the maximum of that, it is guaranteed to be $\geq 0$.
1. Following the below references, we impose priors on quantiles of the GEV distribution. That is, rather than try to use domain knowledge to reason about $\mu,\sigma,\xi$ directly, we instead use domain knowledge to reason about the quantiles of the resulting distribution. We apply LogNormal priors so that a qiven return level is $\geq 0$ and can potentially be very large. The specific LogNormal parameters are chosen based on weak intuition about storm surges along the US Atlantic and Gulf Coast.

> Coles, S. G., & Tawn, J. A. (1996). A Bayesian analysis of extreme rainfall data. Journal of the Royal Statistical Society: Series C (Applied Statistics), 45(4), 463–478. https://doi.org/10.2307/2986068

> Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC. Retrieved from http://ebookcentral.proquest.com/lib/rice/detail.action?docID=4312572
"

# ╔═╡ 6194ae48-258d-4649-9b0a-802393a9e6fa
GEV_priors = NorfolkFloods.get_GEV_priors();

# ╔═╡ 3e525f8f-e94e-4b4a-9ef7-b8ebed4eeb81
md"""
### Visualize Priors

The prior model described above defines a prior distribution over the 10, 100, and 1000 year floods.
This is an alternative to holding prior beliefs over the parameters of a GEV distribution explicitly.
Let's look at this prior!

Formally, this plot is a *prior predictive check* that shows the distribution of three quantities that depend on the parameters $\theta := (\mu, \sigma, \xi)$
"""

# ╔═╡ 14179ad7-74ed-4e31-a6ef-c50c1dd9ba6b
function plot_rl_priors(priors)
    p = plot(
        xlims = (0, 30),
        xlabel = "Storm Surge (ft)",
        ylabel = "Probability Density",
        yticks=[],
		leftmargin=5mm,
    )
    for ((prob, dist), c) in zip(priors, colors)
        return_level = Int(round(1 / (1 - prob)))
        plot!(p, dist, label = "$(return_level) Year Flood", linewidth = 2, color=c)
    end
    return p
end;

# ╔═╡ e2136e56-cbcc-4e4f-9ad9-5437ea3acb01
rl_plot = plot_rl_priors(GEV_priors)

# ╔═╡ d0170cd1-ec50-47b2-b139-b220b36ca2d2
savefig(rl_plot, plotsdir("priors_on_surge_quantiles.pdf"));

# ╔═╡ 528ae93e-bf8e-488d-b1b5-008b4615f4bd
md"""
The above quantiles were reached after some iteration and are by no means *correct*, but they seem subjectively defensible.
See Gelman and Shalizi, 2013, for a philosophical discussion of priors.

> Gelman, A., & Shalizi, C. R. (2013). Philosophy and the practice of Bayesian statistics. British Journal of Mathematical and Statistical Psychology, 66(1), 8–38. https://doi.org/10.1111/j.2044-8317.2011.02037.x

### Build Model

Now let's build our prior model and draw some samples from it
"""

# ╔═╡ dcdd9894-3372-4058-931b-4b6d20500304
prior_model = NorfolkFloods.GEVModel([missing]);

# ╔═╡ 8e7d2aa0-65e8-4441-8a12-9ebd5150041e
md"""
### Sample Model

Next, we want to sample from our prior model in order to check whether the logical implications of our prior beliefs are consistent with what we know about the world.
"""

# ╔═╡ 303563bc-3d38-4589-aedf-2c10b76ce4fb
prior = NorfolkFloods.get_fits(prior_model, "prior_model", n_samples; n_chains = n_chains);

# ╔═╡ 9af4f43f-6225-4692-9d09-c8f61a48e55c
md"
The lowest hanging fruit is to check that the chains have converged, and that they are well mixed.
"

# ╔═╡ 3a3f2707-b206-42e1-abcd-8cda98e10153
summarystats(prior)

# ╔═╡ eb623f59-eea1-42d5-aaaa-331dcb27dba0
plot(prior)

# ╔═╡ f8840050-d813-4680-900a-827772ab860d
md"Both look good thus far"

# ╔═╡ 980550d8-3fa2-4a5b-8619-8b2dc8a3821c
md"""
### Fake Data Experiment

A desirable attribute of a model of storm surges is being able to accurately fit data that comes from a known 'true' distribution.
We can check that with fake data!
"""

# ╔═╡ 3c35e618-eb5b-4dd7-ad41-b3f05f926cd5
begin
    fake_dist = GeneralizedExtremeValue(4, 0.5, 0.15)
    fake_data = DataFrame(CSV.File(datadir("raw", "fake_data.csv")))[!, :fake_data]
    fake_posterior = get_fits(GEVModel(fake_data), "fake_data", n_samples)
    fake_yhat = vcat(
        rand.(
            GeneralizedExtremeValue.(
                fake_posterior[:μ],
                fake_posterior[:σ],
                fake_posterior[:ξ],
            ),
            10,
        )...,
    )
    fake_data_plot = plot(
		xlabel = "Truth",
		ylabel = "Estimated",
        aspect_ratio = :equal,
        legend = :topleft,
    )
    for rt in [2, 5, 10, 25, 50, 100, 250, 500, 1000]
        q = 1 - 1 / rt
        scatter!(
            fake_data_plot,
			[quantile(fake_dist, q)],
			[quantile(fake_yhat, q)],
            label = "",
            color = :gray,
        )
        annotate!(quantile(fake_yhat, q) + 0.1, quantile(fake_dist, q), (rt, 8, :left))
    end
    Plots.abline!(fake_data_plot, 1, 0, label = "1:1 Line")
end

# ╔═╡ e58f7c67-f44c-48e3-8166-d1b2a055e33b
md"""
By looking at this plot, we can see that in this one case, we are able to recover the truth quite nicely with $n_years data points, just like in our observed data.
We note that our prior distribution has introduced a very slight bias, in which the estimates of return levels for short return periods (<50 years) are very slightly over-estimated, and the estimates for higher return periods (>50 years) are slightly under-estimated.
This is to be expected, and one interpretation of Bayesian inference is as a regularizing device that trades off a small bias in exchange for reduced variance in estimates.
In any case, the differences here are quite small, and we can take this as evidence that our estimation procedure is capable of recovering the true distribution from $n_years data points, though this is no absolute guarantee that it will do so with the data set of interest.
"""

# ╔═╡ 29f5e2b9-f384-4535-9463-8524477f085d
savefig(fake_data_plot, plotsdir("surge_fake_data_experiment.pdf"))

# ╔═╡ 09f1655a-8e69-43d5-a0ec-3d8a37eb381e
md"""
Another prior predictive check is to sample directly from *prior predictive* distribution, which is just
$$p(y) = p(y | \theta) p(\theta)$$.
"""

# ╔═╡ 8f20693a-66be-405a-8418-9c2a19119920
begin
	function plot_predictive(p, chains::Chains, N::Int, label::String)
		predictive = NorfolkFloods.sample_predictive(chains, N)
		predictive_flat = vcat(predictive...)[:]
		bins=range(
			quantile(predictive_flat, 0.001),
			quantile(predictive_flat, 0.999),
			length=1000,
		)
		histogram!(
			p,
			predictive_flat,
			bins=bins,
			label=label,
			xlabel="Storm Surge [ft]",
			ylabel="Probability Density",
			yticks=[],
			left_margin=2.5mm,
		)
		return p
	end
	function plot_predictive(chains::Chains, N::Int, label::String)
		p = plot()
		return plot_predictive(p, chains, N, label)
	end
end;

# ╔═╡ f4e8c3e8-8d4e-4b62-9449-b6b3144dea1b
plot_predictive(prior, n_years, "Prior Predictive Distribution")

# ╔═╡ 4997afcf-36a6-456d-ae61-4e2e43e695fe
md"""
We're not saving this plot because we'll add the posterior predictive to it later!

This looks consistent with what we know about storm surge in the Eastern United States in general, before looking at data in Norfolk, VA, which is one way to conceptualize what a prior should do.
In theory our model allows negative values of storm surge: although the shape parameter of the GEV is restricted to be $\xi \geq 0$, which ensures that the distribution is left bounded, the lower bound of a GEV distribution is $\mu - \sigma/\xi$, which can be $<0$.

If we really wanted to get the left tail right, we could go to something like the Blended GEV (Castro-Camilo et al, 2021).
However, here we're concerned with the right tail of the distribution, and so it's less of an issue.

> Castro-Camilo, D., Huser, R., & Rue, H. (2021). Practical strategies for GEV-based regression models for extremes. ArXiv:2106.13110 [Stat]. Retrieved from http://arxiv.org/abs/2106.13110
"""

# ╔═╡ 8a2d4e8a-a50a-48e8-bfd0-e6a2d79bc7bc
md"
## Posterior

Only now that we have run a few basic checks to ensure that our model seems plausible are we ready to actually look at data!
"

# ╔═╡ 12ed5117-a139-4a9e-bde0-a0b5940f9f70
y = ustrip.(u"ft", obs.surge); # can't pass units to the model

# ╔═╡ c26f644e-d6a6-4b12-9f3b-e5d1f70b0c62
md"let's take a quick look at our data; we will develop a more comprehensive visualization later"

# ╔═╡ 6dc4084b-fe26-4f97-b308-cd733ea42413
plot(obs.year, y, label = "", ylabel = "Storm Surge (ft)", marker=".", markersize=4, markercolor=colors[2], color=colors[1])

# ╔═╡ 7fc96550-d53c-46e2-abb0-969d57f85ce7
md"""
### Fit Posterior

Now let's fit our posterior distribution!
The `get_fits` function will automatically cache our results, so it's easier to re-run when we make changes.
"""

# ╔═╡ 4c3d9d19-7e30-428f-9a9d-94e381163725
posterior = NorfolkFloods.get_fits(
	NorfolkFloods.GEVModel(y),
	"posterior",
	n_samples;
	n_chains=4,
);

# ╔═╡ 8ce9fc8f-026d-4370-b98f-7e2decc95afe
md"""
Now we repeat do the usual diagnostics and visualizations of our sampler
"""

# ╔═╡ 94cb1f9e-4bb1-4ec8-9a0a-5a55a5a9ca1d
summarystats(posterior)

# ╔═╡ 5d72be73-ed70-47ad-94c6-1d8ecce9f4fd
plot(posterior)

# ╔═╡ 47485f3d-7d06-4301-8c75-72c53d8b1d1e
md"Once again this looks good; nothing is jumping out as an issue, though again this nis not a guarantee"

# ╔═╡ 970347bf-ace9-499c-91cc-c59d20d63a19
md"""
### Posterior Predictive

Now let's look at our posterior predictive distribution
"""

# ╔═╡ 60cd3263-2999-41eb-83c8-4d90b50bce07
begin
	predictive = plot_predictive(prior, n_years, "Prior Predictive")
	plot_predictive(predictive, posterior, n_years, "Posterior Predictive")
	savefig(predictive, plotsdir("prior_posterior_predictive.pdf"))
	plot(predictive)
end

# ╔═╡ db56826a-32cd-4cfb-b1b3-34552e9e1237
md"""
We can see a notable sharpening from a relatively wide prior to a more narrow posterior.
We can also note

- values below zero are essentially non-existent in the posterior
- the maximum surge recorded ($(maximum(y)) ft) lies towards the upper tail of the posterior predictive distribution
"""

# ╔═╡ 040306ea-6775-40eb-b5da-b00784d62191
md"""
### Posterior Test Statistics

Now that we have a posterior distribution, we can compute some test statistics to see how the distribution of these statistics computed over samples from our posterior predictive distribution compare to the values in our observed data.

For a more nuanced discussion, see

> Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(2), 389–402. https://doi.org/10.1111/rssa.12378
"""

# ╔═╡ ce710bc3-0669-4355-a2ae-e4a8bded7beb
function plot_test_stat(t, yhat; title="", xlabel=:"")
    t_posterior = [t(yi) for yi in yhat]
    observed = t(y)
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
end;

# ╔═╡ 75e9860f-96f1-4c5a-b3b6-832258e158a5
function MannKendall(x::Vector{<:Real})
    return corkendall([1:length(x);], x)
end;

# ╔═╡ 53aa1b9f-0424-4b21-a574-1f61e2922d54
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

# ╔═╡ 42047266-35ee-42ba-bd51-e33a4afc69a1
yhat = NorfolkFloods.sample_predictive(posterior, n_years);

# ╔═╡ 9e7f15ae-a2cf-4139-834f-1c90b2db01c5
test_plots = plot(
    [
		plot_test_stat(t, yhat, title=title, xlabel=xlabel)
		for (title, xlabel, t) in test_stats
	]...,
    size = (1200, 900),
)

# ╔═╡ 3c88fe49-b0eb-407e-aa84-1ade5ede4fc6
savefig(test_plots, plotsdir("surge_test_statistics.pdf"))

# ╔═╡ b72bd68b-e1cc-44e6-9086-03b89329e09d
md"""
We can see that our estimates (histogram and density) are broadly quite consistent with the observed data (vertical green line). The main discrepancies are:

- We seem to consistently underestimate the (absolute value of) partial autocorrelation, which makes sense: we fit an IID model, which is an approximation of the time dependence of the actual climate. There is some evidence for time dependence in this dataset, but the difference is not huge.
- We have a non-zero (albeit very low!) probability of a storm surge in excess of 30 feet. This is not very realistic. However, since we will be taking the convolution of hazard with a damage function that is S-shaped (even if the house is fully submerged, damage cannot be more than the house value), it doesn't matter *for our purposes here*  whether a flood is 30 or 60 feet, so this is not a fatal flaw
"""

# ╔═╡ Cell order:
# ╟─fe66d55d-0565-4740-b0a4-ce82d69cc0e9
# ╠═f0b9ecec-e017-11eb-3f88-c72476818de3
# ╠═38e43be9-d161-42ef-bd06-5649f52fce97
# ╟─e502336a-d13e-468d-8ff1-b0404e20883e
# ╟─b779c4e8-4c80-4958-be3c-28f589555cc3
# ╠═760742a2-da3e-4ce4-a565-4f54c36ec005
# ╟─dd48b4f1-37b5-4698-8b15-bf709ce32dc3
# ╠═56067de6-2d51-48ec-afc3-c42b94642508
# ╠═69b193fa-0d8a-4c4e-81ad-17a776ba7e68
# ╠═0ffeb315-7b5f-4cf4-986e-daae99eecece
# ╠═190eb3a5-91ec-4a40-a3f2-f816daee1e8f
# ╠═c08683e3-718d-400f-9cfc-564d095a3cac
# ╟─6d211cc4-1c17-4788-be8d-972bfe7d3335
# ╠═6194ae48-258d-4649-9b0a-802393a9e6fa
# ╟─3e525f8f-e94e-4b4a-9ef7-b8ebed4eeb81
# ╠═14179ad7-74ed-4e31-a6ef-c50c1dd9ba6b
# ╠═e2136e56-cbcc-4e4f-9ad9-5437ea3acb01
# ╠═d0170cd1-ec50-47b2-b139-b220b36ca2d2
# ╟─528ae93e-bf8e-488d-b1b5-008b4615f4bd
# ╠═dcdd9894-3372-4058-931b-4b6d20500304
# ╟─8e7d2aa0-65e8-4441-8a12-9ebd5150041e
# ╠═303563bc-3d38-4589-aedf-2c10b76ce4fb
# ╟─9af4f43f-6225-4692-9d09-c8f61a48e55c
# ╠═3a3f2707-b206-42e1-abcd-8cda98e10153
# ╠═eb623f59-eea1-42d5-aaaa-331dcb27dba0
# ╟─f8840050-d813-4680-900a-827772ab860d
# ╟─980550d8-3fa2-4a5b-8619-8b2dc8a3821c
# ╠═3c35e618-eb5b-4dd7-ad41-b3f05f926cd5
# ╟─e58f7c67-f44c-48e3-8166-d1b2a055e33b
# ╠═29f5e2b9-f384-4535-9463-8524477f085d
# ╟─09f1655a-8e69-43d5-a0ec-3d8a37eb381e
# ╠═8f20693a-66be-405a-8418-9c2a19119920
# ╠═f4e8c3e8-8d4e-4b62-9449-b6b3144dea1b
# ╟─4997afcf-36a6-456d-ae61-4e2e43e695fe
# ╟─8a2d4e8a-a50a-48e8-bfd0-e6a2d79bc7bc
# ╠═12ed5117-a139-4a9e-bde0-a0b5940f9f70
# ╟─c26f644e-d6a6-4b12-9f3b-e5d1f70b0c62
# ╠═6dc4084b-fe26-4f97-b308-cd733ea42413
# ╟─7fc96550-d53c-46e2-abb0-969d57f85ce7
# ╠═4c3d9d19-7e30-428f-9a9d-94e381163725
# ╟─8ce9fc8f-026d-4370-b98f-7e2decc95afe
# ╠═94cb1f9e-4bb1-4ec8-9a0a-5a55a5a9ca1d
# ╠═5d72be73-ed70-47ad-94c6-1d8ecce9f4fd
# ╟─47485f3d-7d06-4301-8c75-72c53d8b1d1e
# ╟─970347bf-ace9-499c-91cc-c59d20d63a19
# ╠═60cd3263-2999-41eb-83c8-4d90b50bce07
# ╟─db56826a-32cd-4cfb-b1b3-34552e9e1237
# ╟─040306ea-6775-40eb-b5da-b00784d62191
# ╠═ce710bc3-0669-4355-a2ae-e4a8bded7beb
# ╠═75e9860f-96f1-4c5a-b3b6-832258e158a5
# ╠═53aa1b9f-0424-4b21-a574-1f61e2922d54
# ╠═42047266-35ee-42ba-bd51-e33a4afc69a1
# ╠═9e7f15ae-a2cf-4139-834f-1c90b2db01c5
# ╠═3c88fe49-b0eb-407e-aa84-1ade5ede4fc6
# ╟─b72bd68b-e1cc-44e6-9086-03b89329e09d
