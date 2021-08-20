import CSV
import DataFrames
import Dates
import Downloads
import DrWatson
import ProgressBars
import TimeSeries
import Unitful

"data structure to hold hourly data"
struct HourlyGageRecord
    data::TimeSeries.TimeArray
    gage_name::String
    gage_id::String
end

"""
Get the URL for the NOAA file containing tides and currents data

Some things have been hard-coded:
    datum=NAVD
    station=8638610
    time_zone=GMT
    units=metric
"""
function get_noaa_url(year::Int)
    return "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=hourly_height&application=NOS.COOPS.TAC.WL&begin_date=$(year)0101&end_date=$(year)1231&datum=NAVD&station=8638610&time_zone=GMT&units=metric&format=csv"
end

"Get the path of local CSV file where to save a single year of data"
function get_csv_fname(year::Int)
    cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed", "gauge"))
    if !isdir(cache_dir)
        mkdir(cache_dir)
    end
    return joinpath(cache_dir, "$(year).csv")
end

"Parse one of the CSV files downloaded from NOAA"
function parse_noaa_csv_file(fname)
    surge = DataFrames.DataFrame(CSV.File(fname))[[1, 2]]
    DataFrames.rename!(surge, [:datetime, :sl])
    surge[!, :datetime] =
        map(x -> Dates.DateTime(x, "YYYY-mm-dd HH:MM"), surge[!, :datetime])
    surge[!, :sl] = surge[!, :sl]Unitful.u"m" # explicitly specify unit to avoid issues
    surge[!, :sl] = Unitful.uconvert.(Unitful.u"ft", surge[!, :sl])
    DataFrames.dropmissing!(surge)
    return surge
end

"Get a DataFrame of hourly sea level readings"
function get_hourly_obs(year::Int)
    fname = get_csv_fname(year)
    try
        ts = parse_noaa_csv_file(fname)
        return ts
    catch
        url = get_noaa_url(year)
        @info "downloading data for $year"
        Downloads.download(url, fname)
        ts = parse_noaa_csv_file(fname)
        return ts
    end
end

"Get all hourly sea level data"
function produce_norfolk_hourly()
    years = 1928:2020
    df = vcat([get_hourly_obs(year) for year in ProgressBars.tqdm(years)]...)
    hourly = TimeSeries.TimeArray(df, timestamp = :datetime)[:sl]
    return HourlyGageRecord(hourly, "Sewell's Point, VA", "8638610")
end

"convenience function that will save the hourly data"
function get_norfolk_hourly(; overwrite::Bool = false)

    cache_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "data", "processed"))
    cachename = joinpath(cache_dir, "hourly_surge.jld2")

    try
        @assert !overwrite
        hourly = DrWatson.load(cachename, "hourly")
        return hourly
    catch err
        hourly = produce_norfolk_hourly()
        DrWatson.wsave(cachename, Dict("hourly" => hourly))
        return hourly
    end
end
