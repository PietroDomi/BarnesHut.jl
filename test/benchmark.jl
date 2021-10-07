using BarnesHut
using Logging

io = open("test/logs/benchmark_log_dsba_"*string(now())*".log","w+")
logger = SimpleLogger(io)
global_logger(logger)

f = open("test/data/benchamrk_time_"*string(now())*".txt","w")
write(f,"T,N,Alg,θ,time")
N = [10,50,100,500,1000,5000]
T = 24*7

for n in N
    @info "Creating galaxy"
    galaxy = galaxy_start(n,[0.,0.],[0.,0.],10.);
    @info "Simulating a galaxy with N=$n and T=$T\n"

    @info "Starting brute simualation"
    hist_brute = @timed simulationBrute(galaxy,T,3600.);
    write(f,"$T,$n,BF,_,$(hist_brute.time)")
    @info "Brute Simulation in $(hist_brute.time)\n"
    
    @info "Starting tree simulation with θ=0.5"
    hist_tree1 = @timed simulationTree(galaxy,T,3600.,0.4,false);
    write(f,"$T,$n,BH,0.4,$(hist_brute.time)")
    @info "Tree Simulation in $(hist_tree1.time)\n"
    
    @info "Starting tree simulation with θ=0.1"
    hist_tree2 = @timed simulationTree(galaxy,T,3600.,0.8,false);
    write(f,"$T,$n,BH,0.8,$(hist_tree2.time)")
    @info "Tree Simulation in $(hist_tree2.time)\n"

end

@info "End simulation"
flush(io)
close(io)
close(f)