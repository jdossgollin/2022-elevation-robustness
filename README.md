# 2021-elevation-robustness

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project as described in the paper

```bibtex
TBD
```

The code is authored by James Doss-Gollin.

## Setup

To (locally) reproduce this project, do the following:

1. Download this code base.
1. Open a Julia REPL (this code has only been tested with Julia 1.6.4 LTS!) and do the following:

   ```julia
   using Pkg
   Pkg.activate(".")
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
1. Add some text linking to [ABC](https://statmodeling.stat.columbia.edu/2021/11/15/simulation-based-inference-and-approximate-bayesian-computation-in-ecology-and-population-genetics/)
