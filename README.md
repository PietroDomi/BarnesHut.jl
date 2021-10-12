# Barnes Hut Simulation Algorithm
 Implementation of the BarnesHut algorithm in Julia, with benchmark comparisons with the brute force algorithm.

 This is a `Julia` project I created for another course, it is my first attempt with the language, so there's a lot of room for improvement.
 
 To know more about the algorithm: [original paper](https://doi.org/10.1038%2F324446a0); [wikipedia](https://en.wikipedia.org/wiki/Barnes–Hut_simulation); [example](https://jheer.github.io/barnes-hut/).

## Repo structure
You'll find the module and the function under the folder `src/`, while you can use the scripts in `test/` for a quick run.

## Quick start
```
> cd BarnesHut.jl

> julia --project=.

# Go to pkg to download necessary packages
(BarnesHut) pkg> instantiate

julia> include("test/*ANYSCRIPT*.jl")
```