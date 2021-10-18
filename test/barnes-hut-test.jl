using BarnesHut
using Plots

println("Simulating for galaxy...")

galaxy = galaxy_start(10,[0.,0.],[0.,0.],10.)

gal_hist = simulation_tree(galaxy,24,3600.,1.0)

gal_anim = build_animation(gal_hist,[-10.,10.],[-10.,10.];df=2,time_unit="h")

gif(gal_anim, "test-gifs/galaxy-tree-n10-t24-fps60.gif", fps=60)

print("Barnes Hut tested")