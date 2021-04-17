#=
some functions for deailng with units
=#

import Unitful

const meter = Unitful.u"m"
const ft = Unitful.u"ft"

function to_feet(length::T) where {T<:Unitful.Length}
    return Unitful.ustrip(Unitful.uconvert(ft, length))
end
function to_sq_feet(area::T) where {T<:Unitful.Area}
    return Unitful.ustrip(Unitful.uconvert(ft^2, area))
end
