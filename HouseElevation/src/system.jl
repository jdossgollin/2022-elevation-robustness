using JLD2
using Unitful

"""
Define the metrics
"""
struct Outcome{T<:Real}
    upfront_cost_usd::T
    led_usd::T
end

"""
Get all the inputs for the system model
"""
function get_system_model(
    syear::Int,
    eyear::Int;
    fits=get_norfolk_posterior(),
    house_floor_area=1500u"ft^2",
    elevation_init=7u"ft",
    discount_rate=0.02,
    house_value_usd=200_000,
)

    # first, get some functions
    ead_fn = get_expected_annual_damage_emulator(:hazus) # returns cost in $USD as function of MSL
    elevation_cost_fn = get_elevation_cost_function() # cost of elevating as a function of Δh and area

    # next, get some constants
    years = syear:eyear
    N = length(years) # only need to calculate this once

    # weight applied to EAD each year -- compute this once (same for all SOWs)
    Γ = (1 - discount_rate) .^ (years .- syear)

    function f(si::BRICKSimulation, xj::T) where {T<:Unitful.Length}

        # get the house height
        house_elevation = elevation_init + xj # house relative to gauge (zero)

        # outcome 1: up front cost
        upfront_cost_usd = elevation_cost_fn(xj, house_floor_area)

        # NPV damages
        clearance = house_elevation .- si.lsl # house height above MSL as a function of time
        ead_prop = ead_fn.(clearance) # EAD as % of house value, indexed by Time 
        led_prop = sum(ead_prop .* Γ) # lifetime expected damages as proportion of house value
        led_usd = led_prop * house_value_usd # lifetime expected damages as present USD

        # outcome 3: total cost
        total_cost_usd = led_usd + upfront_cost_usd

        # make sure it's the right length and type
        return Outcome(upfront_cost_usd, led_usd)
    end

    # we return **the function**!
    return f
end

"""
Calculate the outcomes for all SOW x Decision pairs

    - `f` is your system model. It should have the form `f(si, c, xj)`.
    - `x` is a discretized decision space
    - `s` is a vector of SOWs
"""
function exhaustive_exploration(
    f, s::Vector{<:BRICKSimulation}, x::Vector{<:Unitful.Length}
)
    I = length(s) # number of states of world
    J = length(x) # number of decisions

    # index metrics by (sow, decision)
    u = hcat([[f(s[i], x[j]) for i in 1:I] for j in 1:J]...)

    return u
end

function get_outcomes(
    x::Vector{<:Unitful.Length};
    syear::Int=2022,
    eyear::Int=2072,
    fits=get_norfolk_posterior(),
    house_floor_area::A=1500u"ft^2",
    elevation_init::L=7u"ft",
    discount_rate::T=0.02,
    house_value_usd::T=200_000.0,
    overwrite::Bool=false,
) where {T<:Real,L<:Unitful.Length,A<:Unitful.Area}
    xmin = minimum(x)
    xmax = maximum(x)
    Nx = length(x)

    # we'll return this
    s = get_norfolk_brick(; syear=syear, eyear=eyear)

    # define a file if we try to load it in
    fname = data_dir(
        "processed",
        (
            "syear_$syear" *
            "eyear_$eyear" *
            "house_floor_area_$house_floor_area" *
            "elevation_init_$elevation_init" *
            "discount_rate_$discount_rate" *
            "house_value_usd_$house_value_usd" *
            "xmin_$xmin" *
            "xmax_$xmax" *
            "Nx_$Nx" *
            ".jld2"
        ),
    )

    # if we to load from file and it works, just do that!
    if !overwrite
        try
            u = JLD2.load(fname, "u")
            return u, s, x # our work is done!
        catch
        end
    end

    # otherwise, need to build it out
    f = get_system_model(
        syear,
        eyear;
        fits=fits,
        house_floor_area=house_floor_area,
        elevation_init=elevation_init,
        discount_rate=discount_rate,
        house_value_usd=house_value_usd,
    )
    u = exhaustive_exploration(f, s, x)

    # save for next time!
    JLD2.save(fname, Dict("u" => u))

    return u, s, x
end
