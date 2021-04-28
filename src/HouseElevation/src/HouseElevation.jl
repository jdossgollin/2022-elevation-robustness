module HouseElevation

include("util.jl")
include("core.jl")
include("deterministic_depth_damage.jl")

export HouseStructure, copy, elevation_cost, depth_damage_frac, depth_damage

end # module
