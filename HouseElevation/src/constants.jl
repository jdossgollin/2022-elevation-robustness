const surge_years = collect(1928:2021) # the years for which we will consider storm surge

const home_dir = dirname(dirname(dirname(@__FILE__))) # the file

"""Get the directory in which data is located"""
datadir(x...) = joinpath(home_dir, "data", x...)

const historic_norfolk_storms = [
    (name="Chesapeake-Potomac", year=1933, is_tc=true),
    (name="Isabel", year=2003, is_tc=true),
    (name="Irene", year=2015, is_tc=true),
    (name="Nor'Ida", year=2009, is_tc=false),
    (name="Ash Wednesday", year=1962, is_tc=false),
]