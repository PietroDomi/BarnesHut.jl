using BarnesHut

println("Starting Galaxy clash simulation")

N = 500
T = 24*30*12
# const plotLimits = nothing
plotLimits = [-20.,40.]
animation = true
fps = 60

galaxy = galaxy_start(N÷2,[0.,0.],[25.,12.],20.)
galaxy2 = galaxy_start(N÷2,[20.,20.],[-25.,-12.],20.)

append!(galaxy,galaxy2)

gal_hist = simulationTree(galaxy,T,3600.,1.)

if animation
    using Plots
    anim_tree = build_animation(gal_hist,plotLimits,plotLimits,24,clash=true,N=N÷2)
    println("Animation computed.\nBuilding gif...")
    gif(anim_tree, "test/gifs/galaxy-clash-n$N-t$T-fps$fps.gif", fps=fps);
end
