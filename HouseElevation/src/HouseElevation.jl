module HouseElevation

include("constants.jl")
include("lsl.jl") # code to read the sea level data
include("surge-bayes.jl")

export data_dir,
    plots_dir,
    TidesAndCurrentsRecord,
    get_annual,
    get_posterior,
    StationaryGEV

end # module
