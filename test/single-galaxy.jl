# if length(LOAD_PATH) == 3
#     # add here the path to the working directory of the module
#     push!(LOAD_PATH,"C:/Users/PietroDomi/OneDrive - Università Commerciale Luigi Bocconi/Documenti/GitHub/barnes-hut-julia")
# end

# using Plots
# using BenchmarkTools
# using BarnesHut

println("Starting Single Galaxy simulation")

N = 500
T = 24*30*4
# const plotLimits = nothing
plotLimits = [-20.,20.]
fps = 60

galaxy = galaxy_start(N÷2,[0.,0.],[0.1,0.1],10.)
# galaxy = cloud_start(N,[0.,0.])


# display(BHNbody.viewStart(galaxy))

gal_hist = simulationTree(galaxy,T,3600.,0.5,false)

anim_tree = build_animation(gal_hist,plotLimits,plotLimits,24);

println("Animation computed.\nBuilding gif...")

gif(anim_tree, "gifs/galaxy-tree-n$N-t$T-fps$fps.gif", fps=fps);
