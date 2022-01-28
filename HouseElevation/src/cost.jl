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
function fit_expected_damage_emulator(key::Symbol=:hazus; N::Int=1_000_000)
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


# constants
elevation_thresholds = [0 5.0 8.5 12.0 14.0][:] # piecewise linear
elevation_rates = [80.36 82.5 86.25 103.75 113.75][:] # cost /ft / ft^2

# interpolation
elevation_itp = LinearInterpolation(elevation_thresholds, elevation_rates)

# user-facing function
function elevation_cost(A::T2, Δh::T1) where {T1<:Unitful.Length, T2<:Unitful.Area}
    area_ft2 = ustrip(u"ft^2", A)
    base_cost = (10000 + 300 + 470 + 4300 + 2175 + 3500) # in USD
    if Δh < 0.0u"ft"
        throw(DomainError(Δh, "Cannot lower the house"))
    elseif Δh ≈ 0.0u"ft"
        cost = 0.0
    elseif 0.0u"ft" < Δh <= 14.0u"ft"
        rate = elevation_itp(ustrip(u"ft", Δh))
        cost = base_cost + area_ft2 * rate
    else
        throw(DomainError(Δh, "Cannot elevate >14ft"))
    end
    return cost
end