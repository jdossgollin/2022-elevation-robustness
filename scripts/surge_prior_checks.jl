#=
Prior checks
=#

using DrWatson
@quickactivate "2021-elevation-robustness"

using Turing, Distributions, DynamicHMC

using Unitful
using Plots, StatsPlots
using NorfolkFloods

priors_on_rl = [
    (0.9, Normal(5, 2)), # 10 year flood
    (0.99, Normal(7, 3)), # 100 year flood
    (0.999, Normal(10, 4)), # 1000 year flood
]

"""
The full GEV model for storm surge:

Assume that the GEV process is stationary. Instead of applying priors over the parameters, use priors over quantiles of the outcome distribution. See Coles & Tawn (1996) and Stephenson (2015). This can be implemented more straightforwardly in Turing.

The specific values of (α, θ) are given above

Coles, S. G., & Tawn, J. A. (1996). A Bayesian analysis of extreme rainfall data. Journal of the Royal Statistical Society: Series C (Applied Statistics), 45(4), 463–478. https://doi.org/10.2307/2986068

Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC. Retrieved from http://ebookcentral.proquest.com/lib/rice/detail.action?docID=4312572
"""
@model function StationaryGEV_quantile_priors(y)

    # weak priors
    μ ~ Flat()
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
prior = sample(prior_model, DynamicNUTS(), 10_000; drop_warmup = true)
plot(prior, [:μ, :σ, :ξ])
yhat_prior = rand.(GeneralizedExtremeValue.(prior[:μ], prior[:σ], prior[:ξ]))[:];
for (prob, prior_dist) in priors_on_rl
    obs = quantile(yhat_prior, prob)
    desired = prior_dist
    @show prob, obs, desired
end

# use historical data for posterior
obs = get_norfolk_annual()
y = ustrip.(u"ft", obs.surge) # observed data, in feet
model = StationaryGEV_quantile_priors(y)
posterior = sample(model, DynamicNUTS(), 10_000, drop_warmup = true)
summarystats(posterior)
plot(posterior)

# posterior predictive distribution
yhat = rand.(GeneralizedExtremeValue.(posterior[:μ], posterior[:σ], posterior[:ξ]))[:];

p1 = plot(obs.year, y, label = "", ylabel = "Storm Surge (ft)")

p2 = histogram(
    yhat,
    label = "Posterior",
    xlabel = "",
    ylabel = "",
    normalize = :pdf,
    orientation = :horizontal,
    linecolor = false,
    fillalpha = 0.5,
)
histogram!(
    p2,
    yhat_prior,
    label = "Prior",
    linecolor = false,
    fillalpha = 0.5,
    orientation = :horizontal,
    normalize = :pdf,
)
plot(p1, p2, link = :y, yrange = (-2.5, 12.5), layout = (1, 2))
