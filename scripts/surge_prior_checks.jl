#=
Prior checks
=#

using DrWatson
@quickactivate "2021-elevation-robustness"

using Turing, Distributions, DynamicHMC

using Unitful
using Plots, StatsPlots, ColorSchemes
using NorfolkFloods

priors_on_rl = [
    (0.9, LogNormal(1.75, 0.5)), # 10 year flood
    (0.99, LogNormal(2, 0.5)), # 100 year flood
    (0.999, LogNormal(2.25, 0.5)), # 1000 year flood
]
p = plot(
    xlims = (0, 30),
    xlabel = "Storm Surge (ft)",
    ylabel = "",
    title = "Prior on Return Levels",
)
for (pr, d) in priors_on_rl
    plot!(p, d, label = "$(Int(round(1 / (1 - pr)))) Year Flood")
end
p
savefig(p, plotsdir("priors_on_surge_quantiles.pdf"))

"""
The full GEV model for storm surge:

Assume that the GEV process is stationary. Instead of applying priors over the parameters, use priors over quantiles of the outcome distribution. See Coles & Tawn (1996) and Stephenson (2015) for related efforts.

Here we constrain ξ >= 0, which requires y to have a lower bound, and  and μ >= 0, which requires that bound to be at least 0. This is reasonable because we know that y >= 0 as the surges are the maxima of residuals, and the residuals have mean zero (by construction; we have subtracted the mean sea level for the year).

Coles, S. G., & Tawn, J. A. (1996). A Bayesian analysis of extreme rainfall data. Journal of the Royal Statistical Society: Series C (Applied Statistics), 45(4), 463–478. https://doi.org/10.2307/2986068

Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC. Retrieved from http://ebookcentral.proquest.com/lib/rice/detail.action?docID=4312572
"""
@model function StationaryGEV_quantile_priors(y)

    # weak priors
    μ ~ FlatPos(0.0)
    σ ~ FlatPos(0.0)
    ξ ~ FlatPos(0.0)

    # we're going to use this distribution repeatedly
    dist = GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for (prob, prior_dist) in priors_on_rl
        rl = quantile(dist, prob)
        Turing.@addlogprob! logpdf(prior_dist, rl)
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end
end

# no data yet for prior
prior_model = StationaryGEV_quantile_priors([missing])

# draw samples from the prior
prior = sample(prior_model, DynamicNUTS(), 50_000; drop_warmup = true)
plot(prior, [:μ, :σ, :ξ])
yhat_prior = rand.(GeneralizedExtremeValue.(prior[:μ], prior[:σ], prior[:ξ]))[:];

# use historical data for posterior
obs = get_norfolk_annual()
y = ustrip.(u"ft", obs.surge) # observed data, in feet
model = StationaryGEV_quantile_priors(y)
posterior = sample(model, DynamicNUTS(), 50_000, drop_warmup = true)
summarystats(posterior)
plot(posterior)

# posterior predictive distribution
yhat = rand.(GeneralizedExtremeValue.(posterior[:μ], posterior[:σ], posterior[:ξ]))[:];

colors = ColorSchemes.tab10

p1 = plot(obs.year, y, label = "", ylabel = "Storm Surge (ft)", title="Observations", marker=".")
p2 = histogram(
    yhat,
    xlabel = "",
    ylabel = "",
    label = "Posterior",
    title="Distribution",
    normalize = :pdf,
    orientation = :horizontal,
    linecolor = false,
    fillalpha = 0.5,
    color = colors[1],
    legend = :bottomright,
)
hline!(p2, [quantile(yhat, 0.99)], color = colors[1], label = "")
histogram!(
    p2,
    yhat_prior,
    label = "Prior",
    linecolor = false,
    fillalpha = 0.5,
    orientation = :horizontal,
    normalize = :pdf,
    color = colors[2],
)
hline!(p2, [quantile(yhat_prior, 0.99)], color = colors[2], label = "")
pp = plot(p1, p2, link = :y, yrange = (-1.25, 10), layout = (1, 2))
pp
savefig(pp, plotsdir("surge_prior_posterior_dists.pdf"))
