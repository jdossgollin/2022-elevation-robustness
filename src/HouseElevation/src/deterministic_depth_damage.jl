import CSV
import DataFrames
import Interpolations
import Unitful

data_dir = abspath(joinpath(@__DIR__, "..", "data"))

"""
Parse the file containing the parameters of the Europa depth-damage model

See [https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R](https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R) for details.
"""
function parse_europa()
    fname = joinpath(data_dir, "fragility_europa.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = dat[!, "depth_feet"]ft
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

"""
Parse the file containing the parameters of the HAZUS depth-damage model

See [https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R](https://github.com/mahkam/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S03_Construction_Cost.R) for details.
"""
function parse_hazus()
    fname = joinpath(data_dir, "fragility_hazus.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = (dat[!, "depth_feet"] .+ 0.0)ft
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

"Generic function to parse the function given a key"
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

"Get the depth-damage interpolation function"
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

"Depth vs damage as a fraction of house value"
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

"Depth vs damage in actual dollars"
function depth_damage(
    house::HouseStructure,
    gauge_depth::T,
    model::Symbol,
) where {T<:Unitful.Length}
    damage_frac = depth_damage_frac(house, gauge_depth, model)
    return damage_frac * house.V
end
