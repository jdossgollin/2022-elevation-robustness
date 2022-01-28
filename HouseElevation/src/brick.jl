using NetCDF
using Unitful

struct BRICKSimulation{I<:Integer,T<:AbstractFloat}
    years::UnitRange{I}
    lsl::Vector{<:Unitful.Length}
    rcp::T
    dynamics::AbstractString
end

function get_norfolk_brick(
    fname::AbstractString=data_dir("raw", "brick", "BRICK.nc");
    syear::Union{Int,Missing}=missing,
    eyear::Union{Int,Missing}=missing,
)

    # read in the file
    rcp_scenarios = NetCDF.ncread(fname, "RCP")
    dynamics = NetCDF.ncread(fname, "dynamics")
    years = Int.(NetCDF.ncread(fname, "time_proj"))
    sims = NetCDF.ncread(fname, "simulation")
    lsl_m = Array{typeof(rcp_scenarios[1]),4}(NetCDF.ncread(fname, "lsl_m"))u"m"
    lsl = uconvert.(u"ft", lsl_m)

    nsim = length(sims)
    nrcp = length(rcp_scenarios)
    ndynam = length(dynamics)

    if ismissing(syear)
        syear = minimum(years)
    end
    idx_syear = findfirst(years .== syear)

    if ismissing(eyear)
        eyear = maximum(years)
    end
    idx_eyear = findfirst(years .== eyear)

    trajectories = [
        BRICKSimulation(
            syear:eyear, lsl[idx_syear:idx_eyear, i, j, k][:], rcp_scenarios[j], dynamics[k]
        ) for i in 1:nsim, j in 1:nrcp, k in 1:ndynam
    ][:]
    return trajectories
end

function get_year_data(t::BRICKSimulation, y::Int)
    idx = findfirst(t.years .== y)
    return t.lsl[idx]
end

function get_year_data(ts::Vector{<:BRICKSimulation}, y::Int)
    idx = findfirst(first(ts).years .== y)
    return [t.lsl[idx] for t in ts]
end