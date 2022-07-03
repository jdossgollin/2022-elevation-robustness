using Distributions
using StatsBase

"""
Invert the weights
"""
function unsort(x_sorted, idx)
    if length(x_sorted) != length(idx)
        throw("lengths must match")
    end
    return [x_sorted[idx[i]] for i in 1:length(idx)]
end

"""
Calculate the weights
"""
function make_weights(prior::T; σ=0.01) where {T<:Distributions.UnivariateDistribution}
    s = get_lsl(; syear=2022, eyear=2100)
    y = ustrip.(u"ft", get_year_data(s, 2100))
    noise = rand(Normal(0, σ), length(y))
    y = y .+ noise
    idx = sortperm(y)
    y_sorted = y[idx]
    cdfs = cdf.(prior, y_sorted)
    w = vcat(first(cdfs), diff(cdfs))
    w = w ./ sum(w)
    w_ordered = unsort(w, idx)
    return StatsBase.Weights(w_ordered)
end
