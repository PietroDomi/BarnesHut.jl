# if length(LOAD_PATH) == 3
#     # add here the path to the working directory of the module
#     push!(LOAD_PATH,"C:/Users/PietroDomi/OneDrive - Università Commerciale Luigi Bocconi/Documenti/GitHub/barnes-hut-julia")
# end

# using Plots
# using BenchmarkTools
using BarnesHut

println("Starting Galaxy clash simulation")

N = 500
T = 24*30*6
# const plotLimits = nothing
plotLimits = [-20.,40.]
fps = 60

galaxy = galaxy_start(N÷2,[0.,0.],[10.,10.],10.)
galaxy2 = galaxy_start(N÷2,[20.,20.],[-10.,-10.],10.)
# galaxy = cloud_start(N,[0.,0.])

append!(galaxy,galaxy2)

# display(BHNbody.viewStart(galaxy))

gal_hist = simulationTree(galaxy,T,3600.,0.5,false)

anim_tree = build_animation(gal_hist,plotLimits,plotLimits,24,clash=true,N=N÷2)

println("Animation computed.\nBuilding gif...")

gif(anim_tree, "gifs/galaxy-clash-n$N-t$T-fps$fps.gif", fps=fps);
