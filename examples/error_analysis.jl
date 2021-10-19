using HDF5
using BarnesHut
using LinearAlgebra
using Statistics
using Plots

N = 1000 .* [i for i in 1:10]
θ = [0.4,0.8,1,1.5,2,5,7,10]

T = 24*7
hist_brute = Dict{String,Array{Float64,3}}()
hist_bh = Dict{String,Array{Float64,3}}()

for n in N
    try
        hist_brute["brute-T$T-N$n"] = h5read("examples/data/output/sim-2021-10-18T19-51-29.399.h5","brute-T$T-N$n")
    catch 
        println(n)
    end
end

for n in N, theta in θ
    try
        hist_bh["barneshut-T$T-N$n-th$theta"] = h5read("examples/data/output/sim-2021-10-18T19-51-29.399.h5","barneshut-T$T-N$n-th$theta")
    catch
        println(n)
    end
end

err = Dict{String,Float64}()
df_err = zeros(80,3)

i = 1

for n in N, theta in θ
    error = 0.0
    println("$n $theta")
    for s in 1:n+1, t in 2:T+1
        error += norm(hist_bh["barneshut-T$T-N$n-th$theta"][t,s,:] - hist_brute["brute-T$T-N$n"][t,s,:])
    end
    err["$n-$theta"] = error
    df_err[i,1] = n
    df_err[i,2] = theta
    df_err[i,3] = error / (n*T)
    i += 1
end

mean_error_theta = zeros(10)

to_plot = zeros(6,10)

plot(to_plot[:,1])

for k in 0:9
    mean_error_theta[k+1] = mean(df_err[k*6+1:k*6+6,3])
    to_plot[:,k+1] = df_err[k*6+1:k*6+6,3]
    plot!(to_plot[:,k+1])
end


open("df_err.csv","w") do f
    print(f, "N,theta,error\n")
    for i in 1:80
        print(f,replace(replace("$(df_err[i,:])\n","["=>""),"]"=>""))
    end
end
