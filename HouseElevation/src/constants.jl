const surge_years = collect(1928:2021) # the years for which we will consider storm surge
const lsl_years = collect(2022:2100)

const home_dir = dirname(dirname(dirname(@__FILE__))) # the file

ensure_dir(dirname) = !isdir(dirname) && mkpath(dirname)

"""Get the directory in which data is located"""
function data_dir(x...)
    dir = joinpath(home_dir, "data", x...)
    ensure_dir(dirname(dir))
    return dir
end

"""Get the directory where we save plots"""
function plots_dir(x...)
    dir = joinpath(home_dir, "plots", x...)
    ensure_dir(dirname(dir))
    return dir
end

"""Clean up all the outputs"""
function clear_cache()
    dirs = [data_dir("processed"), plots_dir()]
    for dir in dirs
        try
            rm(dir; recursive=true)
        catch
        end
    end
end