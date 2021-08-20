module NorfolkFloods

include("hourly.jl")
include("annual.jl")
include("stationary.jl")

export get_norfolk_annual,
    AnnualGageRecord, get_GEV_priors, GEVModel, sample_predictive, get_fits
end # module
