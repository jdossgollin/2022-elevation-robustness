using Base: vcat
using CSV
using DataFrames, DataFramesMeta
using Dates
using Downloads
using ProgressBars
using StatsBase
using Unitful

"""A `WaterLevelStation` stores the basic information about a Tides and Currents station data query"""
Base.@kwdef struct TidesAndCurrentsRecord
    datum::String = "NAVD"
    gage_id::String = "8638610"
    gage_name::String = "Sewells Point, VA"
    units::String = "metric"
    time_zone::String = "GMT"
end

"""
An `HourlyGageRecord` stores all the data from a tide gauge
"""
struct HourlyGageRecord
    t::Vector{<:Dates.DateTime}
    lsl::Vector{<:Unitful.Length}
    stn::TidesAndCurrentsRecord
end

"""
An `AnnualGageRecord` stores mean sea level and annual maximum storm surge for each year
"""
Base.@kwdef struct AnnualGageRecord
    t::Vector{<:Int}
    max_surge::Vector{<:Unitful.Length}
    msl::Vector{<:Unitful.Length}
    stn::TidesAndCurrentsRecord
end;

"""
Get the URL for the NOAA Tides and Currents data; sensible defaults are applied
"""
function get_url(stn::TidesAndCurrentsRecord, year::Int)
    return "https://api.tidesandcurrents.noaa.gov" *
           "/api/prod/datagetter?product=hourly_height&application=NOS.COOPS.TAC.WL" *
           "&begin_date=$(year)0101" *
           "&end_date=$(year)1231" *
           "&datum=$(stn.datum)" *
           "&station=$(stn.gage_id)" *
           "&time_zone=$(stn.time_zone)" *
           "&units=$(stn.units)" *
           "&format=csv"
end

"""
Get the filename in which a given storm surge dataset is to be stored
"""
function get_fname(stn::TidesAndCurrentsRecord, year::Int)
    cache_dir = abspath(datadir("processed", "gage"))
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    return joinpath(
        cache_dir,
        "tidesandcurrents-$(stn.gage_id)-$(year)-$(stn.datum)-$(stn.time_zone)-$(stn.units).csv",
    )
end

"""
Parse a NOAA Tides and Currents file as an `HourlyGageRecord`
"""
function read_tidesandcurrents(stn::TidesAndCurrentsRecord, year::Int)
    fname = get_fname(stn, year)
    date_format = "yyyy-mm-dd HH:MM" # e.g., 1928-01-01 00:00
    surge = @chain CSV.File(fname, dateformat=date_format) begin
        DataFrame
        rename("Date Time" => "datetime", " Water Level" => "lsl")
        select(:datetime, :lsl)
        @transform :lsl = :lsl .* 1u"m"
        @transform :lsl = uconvert.(u"ft", :lsl)
    end
    dropmissing!(surge)
    return HourlyGageRecord(surge[!, :datetime], surge[!, :lsl], stn)
end

"""
Get the observations for an hourly year of data by trying to read from file; if it's not there, download i
"""
function get_hourly(stn::TidesAndCurrentsRecord, year::Int)
    fname = get_fname(stn, year)
    try
        ts = read_tidesandcurrents(stn, year)
        return ts
    catch
        url = get_url(stn, year)
        fdir = dirname(fname)
        !isdir(fdir) && mkdir(fdir)
        Downloads.download(url, fname)
        return read_tidesandcurrents(stn, year)
    end
end

# extend vcat
function Base.vcat(dats::HourlyGageRecord...)
    sdates = [first(dat.t) for dat in dats]
    dats_sorted = dats[sortperm(sdates)]
    t = vcat([dat.t for dat in dats_sorted]...)
    lsl = vcat([dat.lsl for dat in dats_sorted]...)
    stn = first(dats).stn
    return HourlyGageRecord(t, lsl, stn)
end

"""Get the hourly data for all years"""
function get_hourly(stn::TidesAndCurrentsRecord; years=surge_years)
    return vcat([get_hourly(stn, year) for year in years]...)
end

"""Get the annual data for all years"""
function get_annual(stn::TidesAndCurrentsRecord; years=surge_years)
    hourly = [get_hourly(stn, year) for year in years]
    msl = [StatsBase.mean(dat.lsl) for dat in hourly]
    max_lsl = [maximum(dat.lsl) for dat in hourly]
    max_surge = max_lsl .- msl
    return AnnualGageRecord(; t=collect(years), max_surge=max_surge, msl=msl, stn=stn)
end
