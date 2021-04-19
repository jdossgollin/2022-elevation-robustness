#=
This script illustrates how to set a prior on quantiles of the outcome
rather than on the parameters themselves. This is a very useful thing to
do as it is often easier to elicit priors on the outcome variable than on
the parameters themselves.

See the following discussion with Tor Erlend Fjelde (https://github.com/torfjelde)
https://discourse.julialang.org/t/turing-jl-prior-on-quantiles/59510/3
=#

using DrWatson
@quickactivate "2021-elevation-robustness"

using Turing, Distributions, DynamicHMC
using DynamicPPL: @addlogprob!


@model function Demo()

    # define the parameters as vaguely as possible
    μ ~ Turing.Flat()
    σ ~ Turing.FlatPos(0.0)

    # we're going to use this distribution repeatedly
    dist = Normal(μ, σ)

    # define quantiles -- these are deterministic
    q50 = quantile(dist, 0.50)
    q90 = quantile(dist, 0.90)

    # prior distribution on quantiles - pretend this is "from experts"
    DynamicPPL.@addlogprob! logpdf(Normal(3, 0.1), q50) # our prior on the median is about 5
    DynamicPPL.@addlogprob! logpdf(Normal(10, 0.1), q90) # our prior on the 90th percentile is about 10
end

y = zeros(1000)
prior = sample(Demo(), DynamicNUTS(), 100_000)
yhat_prior = rand.(Normal.(prior[:μ], prior[:σ]))[:]
@show quantile(yhat_prior, 0.50) # should be about 3
@show quantile(yhat_prior, 0.90) # should be about 10