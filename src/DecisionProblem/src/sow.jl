"""
Define what a SOW is
"""

using StaticArrays
import StatsBase

struct StateOfWorld{N,T<:Real}
    years::SVector{N,Int}
    msl::SVector{N,T}
    surge::SVector{N,T}
end

function StateOfWorld(;
    years::AbstractVector{Int},
    msl::AbstractVector{T},
    surge::AbstractVector{T},
) where {T<:AbstractFloat}
    length(years) == length(msl) == length(surge) || throw(DimensionMismatch())
    return StateOfWorld{length(years),T}(years, msl, surge)
end
sow1 = StateOfWorld(years = [10, 11, 12], msl = [3.3, 3.4, 3.4], surge = [0.1, -3.1, 5.3])
sow2 = StateOfWorld(years = [5, 11, 3], msl = [3.3, 2.1, 3.4], surge = [0.1, -3.1, 5.3])

"""
A SOWCollection stores a set of SOWs and probabilities corresponding to each
"""
mutable struct SOWCollection{K,N,T<:Real}
    sows::SVector{K,StateOfWorld{N,T}}
    weights::StatsBase.AbstractWeights
end
sc = SOWCollection(SVector(sow1, sow2), StatsBase.uweights(2))
sc.weights = StatsBase.pweights([0.4, 0.6])

x = 1:10
using Statistics: mean
mean(x)
w = StatsBase.pweights(1 ./ x)
mean(x, w)