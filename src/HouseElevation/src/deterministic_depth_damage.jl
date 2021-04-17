#=
Create deterministic depth-damage functions. Could be extended easily enough to stochastic depth-damage functions.
To add another curve:
1. create a function to parse it
2. add an elseif block to parse_data(key)
3. produce a fitted model at the bottom
=#
import CSV
import DataFrames
import Interpolations
import Unitful

data_dir = abspath(joinpath(@__DIR__, "..", "data"))

function parse_europa()
    fname = joinpath(data_dir, "fragility_europa.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = dat[!, "depth_feet"]ft
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

function parse_hazus()
    fname = joinpath(data_dir, "fragility_hazus.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = (dat[!, "depth_feet"] .+ 0.0)ft
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

function parse_data(key)
    if key == :europa
        depth, damage_frac = parse_europa()
    elseif key == :hazus
        depth, damage_frac = parse_hazus()
    else
        throw("Invalid key")
    end
    return depth, damage_frac
end

function get_depth_damage(key)
    depth, damage_frac = parse_data(key)
    prepend!(depth, minimum(depth) - 0.1ft)
    prepend!(damage_frac, 0)
    interp_fn = Interpolations.LinearInterpolation(
        to_feet.(depth),
        damage_frac,
        extrapolation_bc = Interpolations.Flat(),
    )
    damage_fn = function (depth::T) where {T<:Unitful.Length}
        return interp_fn(to_feet(depth))
    end
    return damage_fn
end

# The actual functions that someone will use
depth_damage_frac_europa = get_depth_damage(:europa)
depth_damage_frac_hazus = get_depth_damage(:hazus)

function depth_damage_frac(
    house::HouseStructure,
    gauge_depth::T,
    model::Symbol,
) where {T<:Unitful.Length}
    flood_depth = gauge_depth - house.h
    if model == :europa
        damage_frac = depth_damage_frac_europa(flood_depth)
    elseif model == :hazus
        damage_frac = depth_damage_frac_hazus(flood_depth)
    else
        throw("Invalid model $model")
    end
    return damage_frac
end

function depth_damage(
    house::HouseStructure,
    gauge_depth::T,
    model::Symbol,
) where {T<:Unitful.Length}
    damage_frac = depth_damage_frac(house, gauge_depth, model)
    return damage_frac * house.V
end
