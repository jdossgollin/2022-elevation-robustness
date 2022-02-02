using Unitful

"""
Define the metrics
"""
struct Outcome{T<:Real}
    upfront_cost_usd::T
    led_usd::T
end

function total_cost_usd(u_ij::Outcome)
    return u_ij.upfront_cost_usd + u_ij.led_usd
end

"""
Get all the inputs for the system model
"""
function get_system_model(
    sows::Vector{<:BRICKSimulation};
    surges=get_norfolk_posterior(),
    house_floor_area=1500u"ft^2",
    elevation_init=7u"ft",
    discount_rate=0.02,
    house_value_usd=200_000,
)

    # first, get some functions
    ead_fn = get_expected_annual_damage_emulator(:hazus) # returns cost in $USD as function of MSL
    elevation_cost_fn = get_elevation_cost_function() # cost of elevating as a function of Δh and area

    # next, get some constants
    s1 = first(sows) # the first SOW
    start_year = minimum(s1.years) # we can ID this
    N = length(s1.years) # only need to calculate this once

    # weight applied to EAD each year -- compute this once (same for all SOWs)
    Γ = (1 - discount_rate) .^ (s1.years .- start_year)

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
function exhaustive_exploration(f, s, x)
    I = length(s) # number of states of world
    J = length(x) # number of decisions

    # index metrics by (sow, decision)
    u = hcat([[f(s[i], x[j]) for i in 1:I] for j in 1:J]...)

    return u
end