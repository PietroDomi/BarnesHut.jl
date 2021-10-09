using BarnesHut
using Logging
using Dates

timestamp = replace(string(now()), ":" => "-")

io = open("test/logs/benchmark_log_dsba_"*timestamp*".log","w+")
logger = SimpleLogger(io)
global_logger(logger)

f = open("test/data/benchamrk_time_mean_"*timestamp*".txt","w")
g = open("test/data/benchamrk_time_"*timestamp*".txt","w")
write(f,"T,N,Alg,θ,step_time\n")

N = [10,20,40,80,160,320,640,1280,5120,10240]

T = 24*7

for n in N
    println("\nSimulating for N=$n:\n")
    @info "Creating galaxy"
    galaxy = galaxy_start(n,[0.,0.],[0.,0.],10.);
    @info "Simulating a galaxy with N=$n and T=$T\n"

    @info "Starting brute simualation"
    hist_brute, bench_time_brute = @timed simulationBrute(galaxy,T,3600.;timing=true);
    mean_time_brute = sum(bench_time_brute)/length(bench_time_brute)
    write(g,"Brute, $T, $n\n")
    write(g,"$(bench_time_brute)\n")
    write(f,"$T,$n,BF,_,$(mean_time_brute)\n")
    @info "Avg step Brute Simulation in $(mean_time_brute)\n"
    
    @info "Starting tree simulation with θ=0.4"
    hist_tree1, bench_time_tree1 = simulationTree(galaxy,T,3600.,0.4;timing=true);
    mean_time_tree1 = sum(bench_time_tree1)/length(bench_time_tree1)
    write(f,"$T,$n,BH,0.4,$(mean_time_tree1)\n")
    write(g,"BH, T=$T, n=$n, θ=0.\n4")
    write(g,"$(bench_time_tree1)\n")
    @info "Avg step Tree Simulation in $(mean_time_tree1)\n"
    
    @info "Starting tree simulation with θ=0.8"
    hist_tree2, bench_time_tree2 = @timed simulationTree(galaxy,T,3600.,0.8;timing=true);
    mean_time_tree2 = sum(bench_time_tree2)/length(bench_time_tree2)
    write(f,"$T,$n,BH,0.8,$(mean_time_tree2)\n")
    write(g,"BH, T=$T, n=$n, θ=0.\n8")
    write(g,"$(bench_time_tree2)\n")
    @info "Avg step Tree Simulation in $(mean_time_tree2)\n"

end

@info "End simulation"
flush(io)
close(io)
close(f)
close(g)
