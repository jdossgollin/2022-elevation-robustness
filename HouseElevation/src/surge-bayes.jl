using Distributions
using DynamicHMC
using DynamicPPL
using HDF5
using Latexify
using LaTeXStrings
using MCMCChains
using MCMCChainsStorage
using Optim
using StatsBase
using Turing
using Distributions: InverseGamma

"""
This function converts a Normal distribution to an InverseGamma distribution,
preserving the mean and variance analytically.
"""
function InverseGamma(d::Distributions.Normal)
    α = 2 + d.μ^2 / d.σ^2
    β = (α - 1) * d.μ
    return Distributions.InverseGamma(α, β)
end

# stationary GEV model
gev_priors = [
    (rt=2, dist=InverseGamma(Distributions.Normal(6, 3))),
    (rt=10, dist=InverseGamma(Distributions.Normal(8, 4))),
    (rt=100, dist=InverseGamma(Distributions.Normal(12, 6))),
    (rt=500, dist=InverseGamma(Distributions.Normal(16, 8))),
]

@model function StationaryGEV(y)

    # very weak prior -- need _a_ prior to define variables in Turing
    μ ~ Distributions.TruncatedNormal(0, 25, 0, Inf)
    σ ~ Distributions.TruncatedNormal(0, 10, 0.001, Inf)
    ξ ~ Distributions.TruncatedNormal(0, 2.5, 0, Inf)

    dist = GeneralizedExtremeValue(μ, σ, ξ)

    # implement the prior on quantile levels
    for prior in gev_priors
        xp = 1.0 - 1.0 / prior.rt
        rl = quantile(dist, xp)
        Turing.@addlogprob!(loglikelihood(prior.dist, rl))
    end

    # data model
    if !any(ismissing.(y))
        y .~ dist
    end

    return nothing
end

"""Write a chain to file"""
function write_chain(chain::T, fname::String) where {T<:MCMCChains.Chains}
    mkpath(dirname(fname))
    HDF5.h5open(fname, "w") do f
        write(f, chain)
    end
end

"""Read a chain from file"""
function read_chain(fname::String)
    HDF5.h5open(fname, "r") do f
        read(f, MCMCChains.Chains)
    end
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
        "$(model_name)" * "-$(n_samples)" * "-$(n_chains)" * "-$(drop_warmup)" * ".h5",
    )

    # unless we're overwriting, try to load from file
    if !overwrite
        try
            samples = read_chain(fname)
            return samples
        catch
        end
    end

    # if we're overwriting or reading from file was unsuccessful, sample the model
    sampler = DynamicNUTS() # sampler to use
    x0 = [1.0, 1.0, 1.0] # initial guess
    n_samples = Int(ceil(n_samples / n_chains))

    # get the samples -- no need to thread, it's fast
    samples = mapreduce(
        c -> sample(model, sampler, n_samples; drop_warmup=drop_warmup, init_params=x0),
        chainscat,
        1:n_chains,
    )
    write_chain(samples, fname)
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
    name_replacement = Dict(
        "parameters" => L"\textrm{Parameter}",
        "mean" => L"\textrm{Mean}",
        "std" => L"\textrm{Stdev.}",
        "naive_se" => L"\textrm{Naive SE}",
        "mcse" => L"\textrm{MCSE}",
        "ess" => L"\textrm{ESS}",
        "rhat" => L"$\hat{R}$",
    )
    rename!(df, name_replacement...)
    tex_str = Latexify.latexify(df; env=:table, fmt="%.3f", booktabs=true, index=:subscript)
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
