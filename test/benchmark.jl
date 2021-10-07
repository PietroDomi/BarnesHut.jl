using BarnesHut
using Logging

io = open("test/logs/benchmark_log.log","w+")

logger = SimpleLogger(io)

N = [10,50,100,500,1000,5000,10000]
T = 24*7

global_logger(logger)
for n in N
    @info "Creating galaxy"
    galaxy = galaxy_start(n,[0.,0.],[0.,0.],10.);
    @info "Simulating a galaxy with N=$n and T=$T\n"
    @info "Starting brute simualation"
    hist_brute = @timed simulationBrute(galaxy,T,3600.);
    @info "Brute Simulation in $(hist_brute.time)\n"
    @info "Starting tree simulation with θ=0.5"
    hist_tree1 = @timed simulationTree(galaxy,T,3600.,0.5,false);
    @info "Tree Simulation in $(hist_tree1.time)\n"
    @info "Starting tree simulation with θ=0.1"
    hist_tree2 = @timed simulationTree(galaxy,T,3600.,0.1,false);
    @info "Tree Simulation in $(hist_tree2.time)\n"
end
@info "End simulation"
flush(io)
close(io)
