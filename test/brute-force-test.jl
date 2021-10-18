using BarnesHut
using Plots

println("Simulating for random...")

random = random_start(10,[-10.,10.],[-10.,-10.],5.)

random_hist = simulationBrute(random,24,3600.)

random_anim = build_animation(random_hist,[-50.,50.],[-50.,50.];df=1,time_unit="h")

gif(random_anim, "test-gifs/random-brute-n10-t24-fps60.gif", fps=60)

println("Brute Force tested")