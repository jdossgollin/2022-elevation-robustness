using Distributions
using Turing
using DynamicPPL
using JLD2

data_dir = abspath(joinpath(@__DIR__, "..", "data"))
cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed"))

"""
Get the priors on various return levels
"""
function get_GEV_priors()
    return [
        (0.9, Distributions.LogNormal(1.8, 0.4)), # 10 year flood
        (0.99, Distributions.LogNormal(2, 0.4)), # 100 year flood
        (0.999, Distributions.LogNormal(2.2, 0.4)), # 1000 year flood
    ]
end

gev_priors = get_GEV_priors()

"""
A stationary GEV model. See `get_GEV_priors` for info on the priors.
"""
@model function GEVModel(y)

    # completely flat priors (albeit positive -- see docstring)
    μ ~ Turing.FlatPos(0.0)
    σ ~ Turing.FlatPos(0.0)
    ξ ~ Turing.FlatPos(0.0)

    # we're going to use this distribution repeatedly
    dist = Distributions.GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for (prob, prior_dist) in gev_priors
        rl = Distributions.quantile(dist, prob)
        Turing.@addlogprob! logpdf(prior_dist, rl)
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end
end

"""
Draw samples from the posterior distribution of `fit`.
This will simulate a vector of `N` samples for each  posterior draw in `fit`.
"""
function sample_predictive(fit::Chains, N::Int)
    yhat = [
        rand(Distributions.GeneralizedExtremeValue(μ, σ, ξ), N) for
        (μ, σ, ξ) in zip(fit[:μ], fit[:σ], fit[:ξ])
    ]
    yhat_flat = vcat(yhat...)[:]
    return yhat, yhat_flat
end

"""
No need to re-run computations; this will cache outputs effectively
"""
function get_fits(
    model::DynamicPPL.Model,
    model_name::String,
    n_samples::Int;
    n_chains::Int = 1,
    drop_warmup::Bool = True;
    overwrite::Bool = False,
)
    cachename =
        joinpath(cache_dir, "surge_models", "stationary_$(model_name)_$(n_samples).jld2")
    samples_per_chain = Int(n_samples / n_chains)

    try
        @assert !overwrite
        chains = DrWatson.load(cachename, "chains")
    catch err
        chains = sample(
            model,
            NUTS(),
            MCMCThreads(),
            samples_per_chain,
            n_chains;
            drop_warmup = drop_warmup,
        )
        DrWatson.wsave(cachename, Dict("chains" => chains))
    end
    return chains
end