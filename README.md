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
