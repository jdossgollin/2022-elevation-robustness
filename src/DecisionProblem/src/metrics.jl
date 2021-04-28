using HouseElevation
using DataFrames

"The construction cost is a performance metric"
function construction_cost(x::HouseStructure, l::Lever, sow::StateOfWorld)
    HouseElevation.elevation_cost(x, l.Δh)
end

"NPV damages is a metric"
function npv_damages(
    x::HouseStructure,
    l::Lever,
    sow::StateOfWorld;
    γ = 0.98,
    dmg_model = :europa,
)

    # apply the action: raise the house
    house = copy(x)
    house.h += l.Δh

    # calculate flood damages
    flood = sow.msl .+ sow.surge
    damage = [HouseElevation.depth_damage(house, floodᵢ, dmg_model) for floodᵢ in flood]

    # convert to NPV
    discount = [γ^(Δt) for Δt in (years(sow) .- sow.syear)]
    return sum(discount .* damage)
end

"Compute all performance metrics for a particular decision, on a single `StateOfWorld`"
function evaluate(
    x::HouseStructure,
    l::Lever,
    sow::StateOfWorld;
    γ = 0.98,
    dmg_model = :europa,
)
    metrics = DataFrames.DataFrame(
        :construction_cost => construction_cost(x, l, sow),
        :npv_damages => npv_damages(x, l, sow; γ = γ, dmg_model = dmg_model),
    )
    metrics[!, :total_cost] = metrics[!, :construction_cost] + metrics[!, :npv_damages]
    return metrics
end

"Compute metrics for each SOW"
function evaluate(
    x::HouseStructure,
    l::Lever,
    sc::SOWCollection;
    γ = 0.98,
    dmg_model = :europa,
)
    return vcat([evaluate(x, l, sow; γ = γ, dmg_model = dmg_model) for sow in sc.sows]...)
end

