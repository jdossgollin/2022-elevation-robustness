import CSV
import DataFrames
import DrWatson
import Dates
import JLD2
import TimeSeries
import Unitful
using Statistics: mean

data_dir = abspath(joinpath(@__DIR__, "..", "data"))
cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed"))

norfolk_hourly_fnames = [
    joinpath(data_dir, "1928-01-01_SewellsPoint_1942-08-25.txt"),
    joinpath(data_dir, "1943-09-15_SewellsPoint_2015-12-31.txt"),
]

struct HourlyGageRecord
    data::TimeSeries.TimeArray
    gage_name::String
    gage_id::String
end

"""
Parse a single file containing sea gauge information from NOAA into a TimeArray object
"""
function _parse_noaa_file(fname::Union{String,Vector{String}})

    if typeof(fname) <: Vector
        surge = vcat([CSV.File(f; types = Dict(1 => String)) for f in fname]...)
    else
        surge = CSV.File(fname; types = Dict(1 => String))
    end
    surge = DataFrames.DataFrame(surge)
    surge[!, :datetime] = [
        Dates.Date(date, "YYYYmmdd") + time for
        (date, time) in zip(surge[!, :date], surge[!, :time])
    ]
    surge[!, :sl] = surge[!, :sl]Unitful.u"m" # explicitly specify unit to avoid issues
    return TimeSeries.TimeArray(surge, timestamp = :datetime)[:sl]
end

"""
Read in the Norfolk surge data at Sewell's Point
"""
function produce_norfolk_hourly()
    hourly_surge = _parse_noaa_file(norfolk_hourly_fnames)
    return HourlyGageRecord(hourly_surge, "Sewell's Point", "8638610")
end

"""
Cache the Norfolk surge data if possible
"""
function get_norfolk_hourly(overwrite::Bool = false)
    cachename = joinpath(cache_dir, "hourly_surge.jld2")
    read_raw = true
    try
        hourly = DrWatson.load(cachename, "hourly")
        read_raw = false
    catch err
    end
    if read_raw | overwrite
        hourly = produce_norfolk_hourly()
        DrWatson.wsave(cachename, Dict("hourly" => hourly))
    end
    return hourly
end

"""
Get the return period using a standard formula for the plotting position
"""
function return_period(obs::AbstractVector; a = 0.44)

    @assert length(size(obs)) == 1

    order = sortperm(obs)
    ranks = sortperm(order)
    sample_size = length(obs)

    prob = @. (ranks - a) / (sample_size + 1 - 2 * a)
    return @. 1 / (1 - prob)
end

"""
Create a structure for annual gauge records
"""
struct AnnualGageRecord{I<:Integer,T<:AbstractFloat}
    year::Vector{I}
    surge::Vector{<:Unitful.Length}
    msl::Vector{<:Unitful.Length}
    rt::Vector{T}
    gage_name::String
    gage_id::String

end
function AnnualGageRecord(annual_df::DataFrames.AbstractDataFrame)
    return AnnualGageRecord(
        annual_df[!, :year],
        annual_df[!, :surge],
        annual_df[!, :msl],
        annual_df[!, :rt],
        annual_df[!, :gage_name][1],
        annual_df[!, :gage_id][1],
    )
end
function AnnualGageRecord(hourly::HourlyGageRecord)
    hourly_df = DataFrames.DataFrame(hourly.data)
    hourly_df[!, :year] = [Dates.year(date) for date in hourly_df[!, :timestamp]]
    annual_df = DataFrames.combine(
        DataFrames.groupby(hourly_df, :year),
        :sl => maximum => :surge,
        :sl => mean => :msl,
    )
    annual_df[!, :surge] = annual_df[!, :surge] - annual_df[!, :msl]
    annual_df[!, :rt] = return_period(annual_df[!, :surge])
    annual_df[!, :gage_name] .= hourly.gage_name
    annual_df[!, :gage_id] .= hourly.gage_id
    return AnnualGageRecord(annual_df)
end

"""
It's nice to easily convert the annual gage record to and from a DataFrame
"""
function DataFrame(annual::AnnualGageRecord)
    return DataFrames.DataFrame(
        :year => annual.year,
        :surge => annual.surge,
        :msl => annual.msl,
        :rt => annual.rt,
        :gage_name => annual.gage_name,
        :gage_id => annual.gage_id,
    )
end

"""
Get the annual stage data
"""
function get_norfolk_annual()
    return AnnualGageRecord(get_norfolk_hourly())
end
