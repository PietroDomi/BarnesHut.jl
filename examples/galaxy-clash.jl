using BarnesHut

println("Starting Galaxy clash simulation")

N = 4000
T = 30*12*2
plotLimits = [-60.,60.]
animation = true
cloud = false
fps = 60

if ! cloud
    galaxy = galaxy_start(N÷2,[0.,-20.],[-5000.,1000.],15.)
    galaxy2 = galaxy_start(N÷2,[0.,20.],[5000.,-1000.],15.)
else
    galaxy = cloud_start(N÷2,[-30.,-30.],[20000.,15000.])
    galaxy2 = cloud_start(N÷2,[30.,30.],[-20000.,-15000.])
end

append!(galaxy,galaxy2)

gal_hist = simulation_tree(galaxy,T,24*3600.,4.0)

if animation
    using Plots
    if ! isdir("examples/gifs")
        mkdir("examples/gifs")
    end
    anim_tree = build_animation(gal_hist,plotLimits,plotLimits;df=2,clash=true,N=N÷2,time_unit="d")
    println("Animation computed.\nBuilding gif...")
    gif(anim_tree, "examples/gifs/galaxy-clash-n$N-t$T-fps$fps-2.gif", fps=fps);
end
