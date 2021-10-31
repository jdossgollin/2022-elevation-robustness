# 2021-elevation-robustness

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project as described in the paper

```bibtex
TBD
```

The code is authored by James Doss-Gollin.

## Setup

To (locally) reproduce this project, do the following:

1. Download this code base.
1. Open a Julia console and do the following:

   ```julia
   using Pkg
   Pkg.activate("path/to/this/project")
   Pkg.instantiate()
   include("scripts/main.jl")
   ```

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box.

## Organization

* [`data/`](data/) contains raw data. This data is documented in [`data/raw`](data/raw/).

## Scratchpad (move before publishing)

1. Tony Wong, July 12th:
   1. the `BRICK_NOfastDynamics_SP_20Nov2018.nc` projections are the WK2017 paper’s model projections, redone with Kopp et al 2014 subsidence, extended to 2200, and downscaled to Norfolk instead of NOLA
   1. `BRICK_sewellsPoint_FastDynamics_20Nov2018.nc` is the same thing, but with fast dynamics
1. See original data [on GitHub](https://github.com/scrim-network/local-coastal-flood-risk/tree/master/Data)
