using BarnesHut

println("Starting Single Galaxy simulation")

N = 2000
T = 24*30*12*5
# const plotLimits = nothing
plotLimits = [-20.,20.]
animation = true
fps = 60

galaxy = galaxy_start(N÷2,[0.,0.],[0.0,0.0],10.)
# galaxy = cloud_start(N,[0.,0.],[0.,0.])

gal_hist = simulationTree(galaxy,T,3600.,2.0)

if animation
    using Plots
    anim_tree = build_animation(gal_hist,plotLimits,plotLimits,24);
    println("Animation computed.\nBuilding gif...")
    gif(anim_tree, "test/gifs/galaxy-tree-n$N-t$T-fps$fps.gif", fps=fps);
end