using CSV
using DataFrames
using Interpolations
using Unitful

function parse_europa()
    fname = data_dir("raw", "fragility_europa.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = dat[!, "depth_feet"]u"ft"
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

function parse_hazus()
    fname = data_dir("raw", "fragility_hazus.csv")
    dat = DataFrames.DataFrame(CSV.File(fname))
    depth = (dat[!, "depth_feet"] .+ 0.0)u"ft"
    damage_frac = dat[!, "damage_frac"]
    return depth, damage_frac
end

"""Generic function to get the depth-damage curve given a key `:europa` or `:hazus`"""
function get_depth_damage(key)

    # parse the file
    if key == :europa
        depth, damage_frac = parse_europa()
    elseif key == :hazus
        depth, damage_frac = parse_hazus()
    else
        throw("Invalid key")
    end

    # add a zero to the beginning of the data
    prepend!(depth, minimum(depth) - 0.1u"ft")
    prepend!(damage_frac, 0)

    # interpolate
    interp_fn = Interpolations.LinearInterpolation(
        ustrip.(u"ft", depth), damage_frac; extrapolation_bc=Interpolations.Flat()
    )

    # return *a function*
    damage_fn = function (depth::T) where {T<:Unitful.Length}
        return interp_fn(ustrip.(u"ft", depth))
    end
    return damage_fn
end