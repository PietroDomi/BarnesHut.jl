module BarnesHut

using Plots
using LinearAlgebra
using ProgressBars
using Random

export cloud_start, random_start, galaxy_start, Star
export simulation_brute, simulation_tree
export build_animation

include("utils.jl");
include("quadtree.jl");
include("bruteforce.jl");

println("Building successful")
end # module
