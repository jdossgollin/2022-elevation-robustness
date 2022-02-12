using Distributions
using StatsBase

"""
Calculate the weights
"""
function make_weights(
    s::Vector{<:LSLSim}, prior::T; σ=0.01
) where {T<:Distributions.UnivariateDistribution}
    y = ustrip.(u"ft", get_year_data(s, 2100))
    noise = rand(Normal(0, σ), length(y))
    y = y .+ noise
    idx = sortperm(y)
    y_sorted = y[idx]
    cdfs = cdf.(prior, y_sorted)
    w = vcat(first(cdfs), diff(cdfs))
    w = w ./ sum(w)
    w_ordered = zeros(length(w))

    # get the w in the order that corresponds to the scenarios
    w_ordered = zeros(length(w))
    for i in 1:length(w)
        w_ordered[i] = w[idx[end - i + 1]]
    end
    return StatsBase.Weights(w_ordered)
end