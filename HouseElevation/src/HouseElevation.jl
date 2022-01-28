module HouseElevation

include("constants.jl")
include("lsl.jl") # code to read the sea level data
include("surge-bayes.jl") # code to model
include("brick.jl")

export data_dir,
    plots_dir,
    TidesAndCurrentsRecord,
    get_annual,
    get_posterior,
    StationaryGEV,
    get_norfolk_brick,
    get_year_data

end # module
