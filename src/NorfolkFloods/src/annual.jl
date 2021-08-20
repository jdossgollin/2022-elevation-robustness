using Base: @kwdef
import DataFrames
using Statistics: mean
using DataFrames: DataFrame # to extend



"Create a structure for annual gauge records"
@kwdef struct AnnualGageRecord
    year::Vector{<:Int}
    surge::Vector{<:Unitful.Length}
    msl::Vector{<:Unitful.Length}
    gage_name::String
    gage_id::String
end

function AnnualGageRecord(annual_df::DataFrames.AbstractDataFrame)
    return AnnualGageRecord(
        year = annual_df[!, :year],
        surge = annual_df[!, :surge],
        msl = annual_df[!, :msl],
        gage_name = annual_df[!, :gage_name][1],
        gage_id = annual_df[!, :gage_id][1],
    )
end

function AnnualGageRecord(hourly::HourlyGageRecord)
    hourly_df = DataFrames.DataFrame(hourly.data)
    DataFrames.dropmissing!(hourly_df)
    hourly_df[!, :year] = [Dates.year(date) for date in hourly_df[!, :timestamp]]
    annual_df = DataFrames.combine(
        DataFrames.groupby(hourly_df, :year),
        :sl => maximum => :surge,
        :sl => mean => :msl,
    )
    annual_df[!, :gage_id] .= hourly.gage_id
    annual_df[!, :gage_name] .= hourly.gage_name
    annual_df[!, :surge] = annual_df[!, :surge] - annual_df[!, :msl]
    return AnnualGageRecord(annual_df)
end

"Get the annual stage data"
function get_norfolk_annual()
    return AnnualGageRecord(get_norfolk_hourly())
end