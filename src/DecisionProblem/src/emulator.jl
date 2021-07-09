import Interpolations
import JLD2
import Unitful
import DrWatson
using StatsBase: mean

import HouseElevation
import NorfolkFloods

data_dir = abspath(joinpath(@__DIR__, "..", "data"))
cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed"))

const ft = Unitful.u"ft"

"Get a bunch of synthetic storm surges!"
function draw_surges(n_samples::Int)
    obs = NorfolkFloods.get_norfolk_annual()
    N = Int(ceil(n_samples / 100_000))
    y = Unitful.ustrip.(ft, obs.surge)
    model = NorfolkFloods.GEVModel(y)
    posterior = NorfolkFloods.get_fits(model, "stationary", 100_000; n_chains = 4)
    y_hat = NorfolkFloods.sample_predictive(posterior, N)
    return vcat(y_hat...) .* 1ft
end

"""
Use Monte Carlo approximation to estimate annual expected damages given  house's elevation relative to mean sea level

The expected damages ``d`` as a function of ``h``, the difference between the elevation of the house and mean sea level, is the convolution of the probability density of storm surge ``p(y')`` and the depth-damage function ``d(y' - h)``:

```math
\\mathbb{E}(d | h) = \\int p(y') d(y') dy'
````

This can be approximated as

```math
\\frac{1}{N} \\sum_{i=1}^N d(y^*_i)
````

where ``y^*_1, \\ldots, y^*_N`` are draws from ``p(y')``.
"""
function expected_damage_frac(h::Unitful.Length, s::Vector{<:Unitful.Length}; key = :hazus)
    house = HouseElevation.HouseStructure(2000ft^2, 150_000, h)
    dmg = [HouseElevation.depth_damage_frac(house, si, key) for si in s]
    return mean(dmg)
end

"Fit an interpolation function that estimates expected annual damage given difference between house elevation and MSL"
function fit_expected_damage_emulator(key::Symbol)
    depths = collect(-30:0.5:30)ft
    surge_mc = draw_surges(1_000_000)
    xp_dmg = [expected_damage_frac(depth, surge_mc; key = key) for depth in depths]
    interp_fn = Interpolations.LinearInterpolation(
        HouseElevation.to_feet.(depths),
        xp_dmg,
        extrapolation_bc = Interpolations.Flat(),
    )
    damage_fn = function (h::T) where {T<:Unitful.Length}
        return interp_fn(HouseElevation.to_feet(h))
    end
    return damage_fn
end

"Get an interpolation function that estimates expected annual damage given difference between house elevation and MSL"
function get_expected_damage_emulator(key::Symbol; overwrite::Bool = false)

    valid_keys = [:hazus, :europa]
    @assert(key in valid_keys, "key $key not in $valid_keys")

    cachename = joinpath(cache_dir, "expected_depth_damage_$key.jld2")

    try
        @assert !overwrite
        damage_fn = DrWatson.load(cachename, "damage_fn")
        return damage_fn
    catch err
        damage_fn = fit_expected_damage_emulator(key)
        DrWatson.wsave(cachename, Dict("damage_fn" => damage_fn))
        return damage_fn
    end
end