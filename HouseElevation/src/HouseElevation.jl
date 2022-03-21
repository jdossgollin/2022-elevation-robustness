module HouseElevation

include("constants.jl")
include("lsl.jl") # code to read the sea level data
include("surge-bayes.jl") # code to model
include("io_lsl.jl") # read in the raw data
include("cost.jl")
include("system.jl")
include("aggregate.jl")
include("priors.jl")

export data_dir,
    plots_dir,
    TidesAndCurrentsRecord,
    get_annual_data,
    get_lsl,
    get_year_data,
    get_surge_posterior,
    get_outcomes,
    get_metrics,
    make_weights,
    get_priors

end # module
