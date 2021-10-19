using BarnesHut

println("Starting Galaxy clash simulation")

N = 2000
T = 24*30*12
plotLimits = [-20.,20.]
animation = true
cloud = false
fps = 60

if ! cloud
    galaxy = galaxy_start(N÷2,[0.,-10.],[10000.,0.],10.)
    galaxy2 = galaxy_start(N÷2,[0.,10.],[-10000.,0.],10.)
else
    galaxy = cloud_start(N÷2,[-30.,-30.],[20000.,15000.])
    galaxy2 = cloud_start(N÷2,[30.,30.],[-20000.,-15000.])
end

append!(galaxy,galaxy2)

gal_hist = simulation_tree(galaxy,T,3600.,4.0)

if animation
    using Plots
    anim_tree = build_animation(gal_hist,plotLimits,plotLimits;df=24,clash=true,N=N÷2)
    println("Animation computed.\nBuilding gif...")
    gif(anim_tree, "examples/gifs/galaxy-clash-n$N-t$T-fps$fps.gif", fps=fps);
end
