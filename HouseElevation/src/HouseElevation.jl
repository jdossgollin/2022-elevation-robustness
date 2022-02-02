module HouseElevation

include("constants.jl")
include("lsl.jl") # code to read the sea level data
include("surge-bayes.jl") # code to model
include("brick.jl")
include("cost.jl")
include("system.jl")

export data_dir,
    plots_dir,
    TidesAndCurrentsRecord,
    get_annual_data,
    get_norfolk_brick,
    get_year_data,
    get_norfolk_posterior,
    get_outcomes,
    get_metrics

end # module
