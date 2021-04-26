using Distributions
using Turing
using DynamicPPL

data_dir = abspath(joinpath(@__DIR__, "..", "data"))
cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed"))

import JLD2

"""The GEV model"""
@model function GEVModel(y)

    # completely flat priors (albeit positive -- see docstring)
    μ ~ Turing.FlatPos(0.0)
    σ ~ Turing.FlatPos(0.0)
    ξ ~ Turing.FlatPos(0.0)

    # we're going to use this distribution repeatedly
    dist = Distributions.GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for (prob, prior_dist) in priors_on_rl
        rl = Distributions.quantile(dist, prob)
        Turing.@addlogprob! logpdf(prior_dist, rl)
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end
end

"""Draw samples from the posterior distribution of fit. This will draw one sample for each 
posterior draw in fit. N gives the number of samples per draw, ie the length of the synthetic series.
"""
function sample_predictive(fit::Chains, N::Int)
    yhat = [
        rand(GeneralizedExtremeValue(μ, σ, ξ), N) for
        (μ, σ, ξ) in zip(fit[:μ], fit[:σ], fit[:ξ])
    ]
    yhat_flat = vcat(yhat...)[:]
    return yhat, yhat_flat
end

"""No need to re-run computations; this will cache outputs effectively"""
function get_fit(
    model::DynamicPPL.Model,
    model_name::String,
    n_samples;
    n_chains::Int = 1,
    drop_warmup::Bool = True;
    overwrite::Bool = False,
)
    cachename =
        joinpath(cache_dir, "surge_models", "stationary_$(model_name)_$(n_samples).jld2")
    samples_per_chain = Int(n_samples / n_chains)
    read_raw = true
    try
        chains = DrWatson.load(cachename, "chains")
        read_raw = false
    catch err
    end
    if read_raw | overwrite
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