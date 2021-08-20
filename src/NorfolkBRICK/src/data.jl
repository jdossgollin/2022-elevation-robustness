using NetCDF
using Unitful

data_dir = abspath(joinpath(@__DIR__, "..", "data"))

"""
A single BRICK trajectory has the following information: time (defined by start and end year), local sea level, an RCP scenario, and a BRICK dynamics model (fast or slow)

> Ruckert, K. L., Srikrishnan, V., & Keller, K. (2019). Characterizing the deep uncertainties surrounding coastal flood hazard projections: A case study for Norfolk, VA. Scientific Reports, 9(1), 1–12. https://doi.org/10/ggfsbt =#
"""
struct BRICKtrajectory{I<:Integer,T<:AbstractFloat}
    years::UnitRange{I}
    lsl_m::Vector{<:Unitful.Length}
    rcp::T
    dynamics::AbstractString
end

"""
Parse the naive BRICK sea level simulations. It's very quick so no need to cache. Note that the default value for `fname` assumes the particular directory structure of this project; if you are using this for another project you will need to specify the filename explicitly or else change this argument manually.
"""
function get_norfolk_brick(
    fname::AbstractString = joinpath(data_dir, "BRICK.nc");
    syear::Union{Int,Missing} = missing,
    eyear::Union{Int,Missing} = missing,
)

    @assert isfile(fname) "$fname is not a valid filename"

    rcp_scenarios = NetCDF.ncread(fname, "RCP")
    dynamics = NetCDF.ncread(fname, "dynamics")
    years = Int.(NetCDF.ncread(fname, "time_proj"))
    sims = NetCDF.ncread(fname, "simulation")
    lsl = Array{typeof(rcp_scenarios[1]),4}(NetCDF.ncread(fname, "lsl_m"))Unitful.u"m" # meters

    nsim = length(sims)
    nrcp = length(rcp_scenarios)
    ndynam = length(dynamics)

    if ismissing(syear)
        syear = minimum(years)
    else
        @assert minimum(years) <= syear
    end
    if ismissing(eyear)
        eyear = maximum(years)
    else
        @assert eyear <= maximum(years)
    end

    idx_eyear = findfirst(years .== eyear)
    idx_syear = findfirst(years .== syear)
    trajectories = [
        BRICKtrajectory(
            syear:eyear,
            lsl[idx_syear:idx_eyear, i, j, k][:],
            rcp_scenarios[j],
            dynamics[k],
        ) for i = 1:nsim, j = 1:nrcp, k = 1:ndynam
    ][:]
    return trajectories
end

"Get the data from a particular year"
function get_year_data(t::BRICKtrajectory, y::Int)
    idx = findfirst(t.years .== y)
    return t.lsl_m[idx]
end

"Get the data from a particular year"
function get_year_data(ts::Vector{<:BRICKtrajectory}, y::Int)
    idx = findfirst(first(ts).years .== y)
    return [t.lsl_m[idx] for t in ts]
end