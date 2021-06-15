# if length(LOAD_PATH) == 3
#     # add here the path to the working directory of the module
#     push!(LOAD_PATH,"C:/Users/PietroDomi/OneDrive - Università Commerciale Luigi Bocconi/Documenti/GitHub/barnes-hut-julia")
# end

using Plots
using BarnesHut

Plots.GRBackend()

println("Starting Earth-Sun simulation")

const T = 24*365
# const plotLimits = nothing
const plotLimits = [-20.,20.]
const fps = 60

points = Star[]
BarnesHut.addEarthSun(points,[0.,0.])

es_hist = simulationBrute(points,T,3600.)

es_anim = build_animation(es_hist,plotLimits,plotLimits,24)

println("Animation computed.\nBuilding gif...")

gif(es_anim, "gifs/earthsun-brute-$T-$fps-fps.gif", fps=fps)