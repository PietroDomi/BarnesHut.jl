# if length(LOAD_PATH) == 3
#     # add here the path to the working directory of the module
#     push!(LOAD_PATH,"C:/Users/PietroDomi/OneDrive - Università Commerciale Luigi Bocconi/Documenti/GitHub/barnes-hut-julia")
# end

using Plots
using BarnesHut

# Plots.GRBackend()

println("Starting Earth-Sun simulation")

T = 365*10
N = 20
plotLimits = [-20.,20.]
fps = 60

# points = Star[]
# BarnesHut.addEarthSun(points,[0.,0.])
points = galaxy_start(N,[0.,0.],[0.0,0.0],10.)

es_hist = simulationBrute(points,T,24*3600.)

es_anim = build_animation(es_hist,plotLimits,plotLimits;df=2,time_unit="d");

println("Animation computed.\nBuilding gif...")

gif(es_anim, "examples/gifs/galaxy-brute-n$N-t$T-fps$fps.gif", fps=fps);