module BarnesHut

using Plots
export cloud_start, random_start, galaxy_start, Star
export simulationBrute, simulationTree
export build_animation

include("utils.jl");
include("quadtree.jl");
include("bruteforce.jl");

end # module
