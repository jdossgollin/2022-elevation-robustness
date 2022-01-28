using CSV
using DataFrames
using Interpolations
using JLD2
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

"""Fit the expected damage emulator"""
function fit_expected_damage_emulator(key::Symbol=:hazus; N::Int=10_000)
    clearances = collect(-30:0.5:30)u"ft" # house above MSL
    surge_fits = get_norfolk_posterior()
    surges = sample_predictive_GEV(surge_fits, N)u"ft"

    dmg_fn = get_depth_damage(key)
    xp_dmg = [mean(dmg_fn.(surges .- ci)) for ci in clearances]

    x = ustrip.(u"ft", clearances)
    y = xp_dmg
    interp_fn = Interpolations.LinearInterpolation(
        x, y; extrapolation_bc=Interpolations.Flat()
    )
    damage_fn = function (x::T) where {T<:Unitful.Length}
        return interp_fn(ustrip.(u"ft", x))
    end
    return damage_fn
end

"""Get the expected damage emulator"""
function get_expected_damage_emulator(key::Symbol; overwrite::Bool=false)

    # where to save the emulator
    fname = data_dir("processed", "expected_depth_damage_$key.jld2")

    if !overwrite
        try
            damage_fn = load(fname, "damage_fn")
            return damage_fn
        catch
        end
    end

    damage_fn = fit_expected_damage_emulator(key)
    save(fname, Dict("damage_fn" => damage_fn))
    return damage_fn
end
