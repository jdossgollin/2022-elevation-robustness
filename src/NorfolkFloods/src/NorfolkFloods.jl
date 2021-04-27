module NorfolkFloods

include("surge.jl")
export get_norfolk_hourly,
    get_norfolk_annual,
    HourlyGageRecord,
    AnnualGageRecord,
    DataFrame,
    get_GEV_priors,
    GEVModel,
    sample_predictive,
    get_fits
end # module
