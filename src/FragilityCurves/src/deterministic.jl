using CSV, DataFrames, Interpolations, Unitful

data_dir = abspath(joinpath(@__DIR__, "..", "data"))

"""
Helper function to convert all lengths to feet
"""
function to_feet(length::T) where T <: Unitful.Length
    return Unitful.ustrip(Unitful.uconvert(u"ft", length))
end


function parse_europa()
    fname = joinpath(data_dir, "fragility_europa.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = dat[!, "depth_feet"]u"ft"
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

function parse_hazus()
    fname = joinpath(data_dir, "fragility_hazus.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = (dat[!, "depth_feet"] .+ 0.0)u"ft"
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
    prepend!(depth, minimum(depth) - 1u"mm")
    prepend!(damage_frac, 0)
    interp_fn = Interpolations.LinearInterpolation(
        to_feet.(depth),
        damage_frac,
        extrapolation_bc=Interpolations.Flat(),
    )
    damage_fn = function (depth::T) where T <: Unitful.Length
        return interp_fn(to_feet(depth))
    end
    return damage_fn
end
