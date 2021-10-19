# BarnesHut.jl

 Implementation of the BarnesHut algorithm in Julia for solving the N-body problem, along with benchmark comparisons with the brute force algorithm.

 This is a `Julia` project I created for the course [20602 - Computer Science (Algorithms)](https://didattica.unibocconi.eu/ts/tsn_anteprima.php?cod_ins=20602&anno=2021&IdPag=6164) at UniBocconi, under prof. Feinauer and prof. Pittorino.
 
 To know more about the algorithm: 
 - [original paper](https://doi.org/10.1038%2F324446a0)
 - [wikipedia](https://en.wikipedia.org/wiki/Barnes–Hut_simulation)
 - [example](https://jheer.github.io/barnes-hut/)

## Repo structure
The main module and functions are under the folder `src/`.

There are a bunch of scripts ready to run in `examples/`.

## Quick start
Clone the repository first, then access the folder:
```
> cd BarnesHut.jl

> julia --project=.

# Go to pkg to download necessary packages
(BarnesHut) pkg> instantiate

# (Optional) Test the package 
(BarnesHut) pkg> test
```
Next, you can run any of the examples by typing in the REPL

```
julia> include("examples/YOUREXAMPLE.jl")
```

## Credits
Pietro Dominietto (pietro.dominietto@studbocconi.it)
