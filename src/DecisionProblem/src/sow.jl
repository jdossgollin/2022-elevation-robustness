using Distributions
using Unitful

"""
A StateOfWorld contains all the information needed to evaluate metrics of a particular SOW given a decision
    
The `syear` and `eyear` define the first and last time steps. `msl` and `surge` are vectors of lengths.
"""
struct StateOfWorld
    syear::Int
    eyear::Int
    msl::Vector{<:Unitful.Length}
    surge::Vector{Vector{<:Unitful.Length}}
    StateOfWorld(syear, eyear, msl, surge) =
        length(surge) == length(msl) == (eyear - syear + 1) ?
        new(syear, eyear, msl, surge) : error("not all lengths are equal")
end

"You can pass in a Vector or Tuple of years if you prefer"
function StateOfWorld(years, msl::Vector{<:Unitful.Length}, surge::Vector{<:Unitful.Length})
    syear = Int(minimum(years))
    eyear = Int(maximum(years))
    return StateOfWorld(syear, eyear, msl, surge)
end

"Get the years of a particular SOW"
function years(sow::StateOfWorld)
    return Vector(sow.syear:sow.eyear)
end

"Get the number of years of a particular SOW"
function n_years(sow::StateOfWorld)
    return sow.eyear - sow.syear + 1
end
