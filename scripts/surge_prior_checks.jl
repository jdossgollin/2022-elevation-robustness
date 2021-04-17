#=
Prior checks
=#

using Turing, Distributions, Unitful, Plots, StatsPlots, DynamicHMC, Optim
using NorfolkFloods

"""
Estimate the parameters of a Gamma distribution from quantiles

It is often desirable to elicit information on quantiles rather than moments (mean and variance).
The easiest way to do this is through optimization.

q ∈ [0, 1] gives the quantiles of the distribution
d gives the corresponding values

thus quantile.(Gamma(α, θ), q) ≈ d
"""
function gamma_from_quantile(q::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    @assert length(q) == length(d)
    # define the objective for optim
    function objective(x)
        α, θ = x
        dist = Gamma(α, θ)
        return sum((quantile.(dist, q) .- d) .^ 2)
    end
    lower = zeros(2)
    upper = [Inf, Inf]
    initial_x = [2.0, 2.0]
    inner_optimizer = LBFGS()
    res = optimize(objective, lower, upper, initial_x, Fminbox(inner_optimizer))
    return Optim.minimizer(res)
end

"""
The Jacobian for the prior

See the unnumbered equation following equation 13.7 of

Stephenson, A. (2015). Bayesian inference for extreme value modeling. In D. K. Dey & J. Yan, Extreme value modeling and risk analysis: methods and applications. Philadelphia, PA: CRC Press LLC.
"""
function jacobian_coles(θ, x)
    μ, σ, ξ = θ
    return (σ / ξ^2) * abs(
        sum([
            (-1)^(i + j) * (x[i] * x[j])^(-ξ) * log(x[j] / x[i]) for j = 1:length(x) for
            i = 1:max(j - 1, 1)
        ]),
    )
end

function get_prior(θ, x, α, β, q̃ₚ)
    return jacobian_coles(θ, x) *
           prod([q̃ₚ[i]^(α[i] - 1 * exp(-q̃ₚ[i] / β[i])) for i = 1:3])
end

q = [0.5, 0.9]
d = [59, 72]
α, θ = gamma_from_quantile(q, d)
β = 1 / θ
@show α, β
quantile.(Gamma(α, θ), q) ≈ d

"""
A Bayesian stationary GEV model
"""
@model StationaryGEV(y) = begin

    p = [0.1, 0.01, 0.001]
    x = @. -log(1 - p)

    # priors
    #μ ~ Normal(4, 2)
    #σ ~ Truncated(Normal(0, 1), 0, Inf)
    #ξ ~ Normal(0.1, 0.1)

    # data model
    y .~ GeneralizedExtremeValue(μ, σ, ξ)
end

obs = get_norfolk_annual()
y = ustrip.(u"ft", obs.surge)
model = StationaryGEV(y)
prior = sample(model, Prior(), 5_000)
summarystats(prior)

yhat_prior = rand.(GeneralizedExtremeValue.(prior[:μ], prior[:σ], prior[:ξ]))[:]
p0 = histogram(yhat_prior, orientation = :horizontal, title = "Prior", label = "")

posterior = sample(model, DynamicNUTS(), MCMCThreads(), 5000, 4, drop_warmup = true)
plot(posterior)
yhat = rand.(GeneralizedExtremeValue.(posterior[:μ], posterior[:σ], posterior[:ξ]))[:]

p1 = plot(obs.year, y, label = "", ylabel = "Storm Surge (ft)", title = "Observed")
p2 = histogram(yhat, orientation = :horizontal, title = "Posterior", label = "")
plot(p1, p0, p2, link = :y, yrange = (-2.5, 12.5), layout = (1, 3))

quantile(yhat_prior, 0.001)
quantile(yhat_prior, 0.9)
quantile(yhat_prior, 0.99)
quantile(yhat_prior, 0.999)

quantile(yhat, 0.001)
quantile(yhat, 0.9)
quantile(yhat, 0.99)
quantile(yhat, 0.999)
