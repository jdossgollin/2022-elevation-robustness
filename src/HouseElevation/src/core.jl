#=
Define what we mean by a house
=#

import Interpolations, Unitful

"""
A house structure contains all information needed to calculate costs of elevation and floods
"""
mutable struct HouseStructure{P<:Unitful.Area,T<:Real,L<:Unitful.Length}
    A::P # the house area, in square feet
    V::T # the house value in 2020 $USD 
    h::L # the height of the house relative to the gauge
end


elevation_thresholds = [0 5.0 8.5 12.0 14.0][:] # cost is piecewise linear
elevation_rates = [80.36 82.5 86.25 103.75 113.75][:] # cost per unit increase per area
elevation_itp = Interpolations.LinearInterpolation(elevation_thresholds, elevation_rates)
"""
Get the elevation cost that goes with a particular house folliwing Zarekarizi et al (2020)
This is in turn based on the CLARA model (described in Zarekarizi et al).

Key thing to know is that the cost is piecewise linear and the valid domain of
the input heightening is [0, 14ft]
"""
function elevation_cost(house::HouseStructure, Δh::T) where {T<:Unitful.Length}
    area_ft2 = to_sq_feet(house.A)
    base_cost = (10000 + 300 + 470 + 4300 + 2175 + 3500) # all costs are in dollars
    if Δh < 0.0ft
        throw(DomainError(Δh, "Cannot lower the house"))
    elseif Δh ≈ 0.0ft
        cost = 0.0
    elseif 0.0ft < Δh <= 14.0ft
        rate = elevation_itp(to_feet(Δh))
        cost = base_cost + area_ft2 * rate
    else
        throw(DomainError(Δh, "Cannot elevate >14ft"))
    end
    return cost
end
