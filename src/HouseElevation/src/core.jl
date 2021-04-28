import Base: copy
import Interpolations, Unitful

"All information needed to compute costs for a particular house"
mutable struct HouseStructure{P<:Unitful.Area,T<:Real,L<:Unitful.Length}
    A::P # the house area, in square feet
    V::T # the structural value of the house, in USD2020
    h::L # the height of the house relative to the gauge
end

"copy a house structure"
function copy(hs::HouseStructure)
    return HouseStructure(hs.A, hs.V, hs.h)
end


elevation_thresholds = [0 5.0 8.5 12.0 14.0][:] # cost is piecewise linear
elevation_rates = [80.36 82.5 86.25 103.75 113.75][:] # cost per unit increase per area
elevation_itp = Interpolations.LinearInterpolation(elevation_thresholds, elevation_rates)

"""
    elevation_cost(house, Δh)

Compute the cost of elevating a particular HouseStructure `house` by `Δh`, where `house` is a `HouseStructure` and `Δh` is a length.

This cost function follows Zarekarizi et al (2020), which is in turn based on the CLARA model. The valid domain of `Δh` is [0, 14ft].

> Zarekarizi, M., Srikrishnan, V., & Keller, K. (2020). Neglecting uncertainties biases house-elevation decisions to manage riverine flood risks. Nature Communications, 11(1), 5361. https://doi.org/10.1038/s41467-020-19188-9
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
