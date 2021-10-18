using BarnesHut

println("Starting Single Galaxy simulation")

N = 2000
T = 30#*12*10
# const plotLimits = nothing
plotLimits = [-40.,40.]
animation = true
fps = 60

galaxy = galaxy_start(N,[0.,0.],[0.0,0.0],20.)
# galaxy = cloud_start(N,[0.,0.],[0.,0.])

gal_hist = simulation_tree(galaxy,T,24*3600.,2.0)

if animation
    using Plots
    anim_tree = build_animation(gal_hist,plotLimits,plotLimits);
    println("Animation computed.\nBuilding gif...")
    gif(anim_tree, "test/gifs/galaxy-tree-n$N-t$T-fps$fps.gif", fps=fps÷2);
end