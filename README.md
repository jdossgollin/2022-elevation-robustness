# 2021-elevation-robustness

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
> 2021-elevation-robustness

It is authored by James Doss-Gollin.

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the git-history and may need to be downloaded independently.
1. Open a Julia console and do:

   ```julia
   using Pkg
   Pkg.activate("path/to/this/project")
   Pkg.instantiate()
   include("scripts/main.jl")
   ```

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box.

## Scratchpad

To add the local packages

```julia
Pkg.develop(PackageSpec(path="directory"))
```
