using BarnesHut

println("Starting Earth-Sun simulation")

T = 365*10
N = 20
plotLimits = [-20.,20.]
animation = true
fps = 60

points = Star[]
BarnesHut.add_earth_sun(points,[0.,0.])

es_hist = simulation_brute(points,T,24*3600.)

if animation
    using Plots
    if ! isdir("examples/gifs")
        mkdir("examples/gifs")
    end
    es_anim = build_animation(es_hist,plotLimits,plotLimits;df=2,time_unit="d")
    println("Animation computed.\nBuilding gif...")
    gif(es_anim, "examples/gifs/galaxy-brute-n$N-t$T-fps$fps.gif", fps=fps)
end