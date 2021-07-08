using HouseElevation
using DataFrames
using StatsBase

"Check whether all elements of `x` are equal"
allsame(x) = all(y -> y == first(x), x)

"""
Evaluate a fixed decision on a single state of the world

`x` gives the initial house to be elevated
`l` gives the elver to be evaluated
`sow` describes the state of the world
`γ` gives the discount rate
`dmg_model` gives the depth_damage model used

Returns a `Dict` of metrics
"""
function evaluate(
    x::HouseStructure,
    l::Lever,
    sow::StateOfWorld,
    γ::Real,
    dmg_model::Symbol,
)

    # apply the action: raise the house
    house = copy(x)
    house.h += l.Δh
    construction_cost = elevation_cost(x, l.Δh)

    # calculate flood damages
    flood = sow.msl .+ sow.surge
    damage = [depth_damage(house, fᵢ, dmg_model) for fᵢ in flood]

    # convert to NPV
    discount = [γ^(Δt) for Δt in (years(sow) .- sow.syear)]
    npv_damages = sum(discount .* damage)

    # return a dict of metrics
    return Dict(
        :construction_cost => construction_cost,
        :npv_damages => npv_damages,
        :total_cost => construction_cost + npv_damages,
    )
end

"""
Evaluate a fixed decision on a set of states of the world

Returns a `DataFrame` of metrics indexed by their State of World
"""
function evaluate(
    x::HouseStructure,
    l::Lever,
    sows::Vector{StateOfWorld};
    γ::Real = 0.98,
    dmg_model::Symbol = :europa,
)
    return vcat([DataFrame(evaluate(x, l, sow, γ, dmg_model)) for sow in sows]...)
end

"""
Get the expected performance of a lever, given a set of weights

```julia
expected_performance(metrics, w)
```

where `metrics` is the output from `evaluate()` and `w` are weights
"""
function expected_performance(metrics::DataFrame, w::AbstractWeights)
    @assert nrow(metrics) == length(w)
    mapcols(x -> mean(x, w), metrics)
end
