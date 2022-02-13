using Distributions
using DynamicPPL
using Latexify
using LaTeXStrings
using MCMCChains
using Optim
using StatsBase
using Turing

# stationary GEV model
gev_priors = [
    (rt=2, dist=Distributions.TruncatedNormal(4, 1.5, 0, Inf)),
    (rt=10, dist=Distributions.TruncatedNormal(6, 1.75, 0, Inf)),
    (rt=100, dist=Distributions.TruncatedNormal(10, 2.25, 0, Inf)),
    (rt=500, dist=Distributions.TruncatedNormal(15, 2.75, 0, Inf)),
]

@model function StationaryGEV(y)

    # have to put a prior to define parameters in Turing
    ξ ~ Distributions.TruncatedNormal(0, 0.5, 0, Inf) # surge has lower bound => ξ>0
    μ ~ Distributions.TruncatedNormal(0, 10, 0, Inf) # lower bound should not be negative
    σ ~ Distributions.TruncatedNormal(0, 4, 0, Inf) # σ > 0 by definition

    dist = Distributions.GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for prior in gev_priors
        xp = 1.0 - 1.0 / prior.rt
        rl = quantile(dist, xp)
        #Turing.@addlogprob!(logpdf(prior.dist, rl))
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end

    return nothing
end

"""
Get the posterior fitted model. Uses caching to avoid resampling when feasible.
"""
function get_posterior(
    model::DynamicPPL.Model,
    model_name::String,
    n_samples::Int;
    n_chains::Int=1,
    drop_warmup::Bool=true,
    overwrite::Bool=false,
)
    fname = data_dir(
        "processed",
        "surge_models",
        "$(model_name)" * "-$(n_samples)" * "-$(n_chains)" * "-$(drop_warmup)" * ".jls",
    )
    samples_per_chain = Int(ceil(n_samples / n_chains))

    # unless we're overwriting, try to load from file
    if !overwrite
        try
            samples = read(fname, MCMCChains.Chains)
            return samples
        catch
        end
    end

    # if we're overwriting or reading from file was unsuccessful, sample the model
    samples = sample(
        model, NUTS(), MCMCThreads(), samples_per_chain, n_chains; drop_warmup=drop_warmup
    )
    mkpath(dirname(fname))
    write(fname, samples)
    return samples
end

"""The Weibull plotting position"""
function weibull_plot_pos(y)
    N = length(y)
    ys = sort(y; rev=false) # sorted values of y
    nxp = xp = [r / (N + 1) for r in 1:N] # exceedance probability
    xp = 1 .- nxp
    return xp, ys
end

"""Save your diagnostics to file"""
function write_diagnostics(fit, fname::String)
    df = DataFrame(MCMCChains.summarize(fit))
    rename!(
        df,
        "parameters" => L"\textrm{Parameters}",
        "mean" => L"\textrm{Mean}",
        "std" => L"\textrm{Stdev.}",
        "naive_se" => L"\textrm{Naive SE}",
        "mcse" => L"\textrm{MCSE}",
        "ess" => L"\textrm{ESS}",
        "rhat" => L"$\hat{R}$",
        "ess_per_sec" => L"\textrm{ESS per second}",
    )
    tex_str = latexify(df; env=:table, fmt="%.2f")
    open(fname, "w") do io
        write(io, tex_str)
    end
    return true
end

"""Sample from the posterior predictive of the GEV"""
function sample_predictive_GEV(fit, N)
    μ = vec(fit[:μ])
    σ = vec(fit[:σ])
    ξ = vec(fit[:ξ])
    idx = sample(1:length(μ), N)
    return rand.(GeneralizedExtremeValue.(μ[idx], σ[idx], ξ[idx]))
end

"""Get the posterior for Norfolk"""
function get_surge_posterior()
    stn = TidesAndCurrentsRecord()
    annual = HouseElevation.get_annual_data(stn)
    surge_ft = Unitful.ustrip.(u"ft", annual.max_surge) # scalarize in ft

    model = StationaryGEV(surge_ft)

    fits = get_posterior(
        model, "surge_posterior", 10_000; n_chains=4, overwrite=false, drop_warmup=true
    )
    return fits
end