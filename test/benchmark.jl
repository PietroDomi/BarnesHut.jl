using BarnesHut
using Logging
using Dates
using HDF5

timestamp = replace(string(now()), ":" => "-")

io = open("test/logs/benchmark_log_dsba_"*timestamp*".log","w+")
logger = SimpleLogger(io)
global_logger(logger)

f = open("test/data/time/benchmark_time_mean_$timestamp.txt","w")
g = open("test/data/time/benchmark_time_$timestamp.txt","w")
write(f,"T,N,Alg,θ,step_time\n")

N = [10,20,40,80,160,320,640,1280,5120,10240]
θ = [0.4,0.8,1,1.5,2,5]

T = 24*7

for n in N
    println("\nSimulating for N=$n:\n")
    @info "Creating galaxy"
    galaxy = galaxy_start(n,[0.,0.],[0.,0.],10.);
    @info "Simulating a galaxy with N=$n and T=$T\n"

    @info "Starting brute simualation"
    hist_brute, bench_time_brute = simulationBrute(galaxy,T,3600.;timing=true);
    # println(typeof(hist_brute))
    mean_time_brute = sum(bench_time_brute)/length(bench_time_brute)
    write(g,"Brute, $T, $n\n")
    write(g,"$(bench_time_brute)\n")
    write(f,"$T,$n,BF,_,$(mean_time_brute)\n")
    @info "Avg step Brute Simulation in $(mean_time_brute)\n"
    h5write("test/data/output/sim-$timestamp.h5","brute-T$T-N$n",hist_brute)
    
    for theta in θ
        @info "Starting tree simulation with θ=$theta"
        hist_tree, bench_time_tree = simulationTree(galaxy,T,3600.,theta;timing=true);
        mean_time_tree = sum(bench_time_tree)/length(bench_time_tree)
        write(f,"$T,$n,BH,$theta,$(mean_time_tree)\n")
        write(g,"BH, T=$T, n=$n, θ=$theta\n")
        write(g,"$(bench_time_tree)\n")
        @info "Avg step Tree Simulation in $(mean_time_tree)\n"
        h5write("test/data/output/sim-$timestamp.h5","barneshut-T$T-N$n-th$theta",hist_tree)
    end
end

@info "End simulation"
flush(io)
close(io)
close(f)
close(g)
