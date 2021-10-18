using BarnesHut
using Logging
using Dates
using HDF5

timestamp = replace(string(now()), ":" => "-")

# io = open("examples/logs/benchmark_log_dsba_$timestamp.log","w+")
# logger = SimpleLogger(io)
# global_logger(logger)

# f = open("examples/data/time/benchmark_time_mean_$timestamp.txt","w")
# g = open("examples/data/time/benchmark_time_$timestamp.txt","w")
# write(f,"T,N,Alg,θ,step_time\n")

# N = 1000 .* [i for i in 1:10]
N = [1000]
θ = [0.4]#,0.8,1,1.5,2,5,7,10]
T = 24*7

open("examples/data/time/benchmark_time_$timestamp.txt","w") do g
    open("examples/data/time/benchmark_time_mean_$timestamp.txt","w") do f
        open("examples/logs/benchmark_log_dsba_$timestamp.log","w+") do io
            logger = SimpleLogger(io)
            global_logger(logger)
            print(f,"T,N,Alg,θ,step_time\n")

            for n in N
                println("\nSimulating for N=$n:\n")
                @info "Creating galaxy"
                galaxy = galaxy_start(n,[0.,0.],[0.,0.],10.);
                @info "Simulating a galaxy with N=$n and T=$T\n"
            
                @info "Starting brute simualation"
                hist_brute, bench_time_brute = simulationBrute(galaxy,T,3600.;timing=true);
                mean_time_brute = sum(bench_time_brute)/length(bench_time_brute)
                print(g,"Brute, $T, $n\n")
                print(g,"$(bench_time_brute)\n")
                print(f,"$T,$n,BF,0.0,$(mean_time_brute)\n")
                @info "Avg step Brute Simulation in $(mean_time_brute)\n"
                h5write("examples/data/output/sim-$timestamp.h5","brute-T$T-N$n",hist_brute)
                
                for theta in θ
                    @info "Starting tree simulation with θ=$theta"
                    hist_tree, bench_time_tree = simulationTree(galaxy,T,3600.,theta;timing=true);
                    mean_time_tree = sum(bench_time_tree)/length(bench_time_tree)
                    print(f,"$T,$n,BH,$theta,$(mean_time_tree)\n")
                    print(g,"BH, T=$T, n=$n, θ=$theta\n")
                    print(g,"$(bench_time_tree)\n")
                    @info "Avg step Tree Simulation in $(mean_time_tree)\n"
                    h5write("examples/data/output/sim-$timestamp.h5","barneshut-T$T-N$n-th$theta",hist_tree)
                end
            end
            @info "End simulation"

        end
    end
end


# flush(io)
# close(io)
# close(f)
# close(g)