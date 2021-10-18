import Base.-
using LinearAlgebra
using ProgressBars
using Random

import Plots: scatter, scatter!

mutable struct Star
    s::Array{Float64,1} # * 10^10
    v::Array{Float64,1}
    m::Float64 # times 10^20
    Star(s,v,m) = new(s,v,m)
end

Star() = Star([0.,0.],[0.,0.],0.)

distance2(a::Star,b::Star) = (a.s[1]-b.s[1])^2+(a.s[2]-b.s[2])^2

function distance(a::Array{Float64,1},b::Array{Float64,1})
    d = 0.
    for i in 1:length(a)
        d += (a[i]-b[i])^2
    end
    return sqrt(d)
end

-(a::Star,b::Star) = [a.s[1]-b.s[1],a.s[2]-b.s[2]]

function newton(a::Star,b::Star;max::Float64=1e18)
    G = 6.674
    r2 = distance2(a,b)
    F = min(G * a.m * b.m / r2, max)  # * 10^9
    return F 
end

function get_cos_sin(a::Array)
    cos = a[1]/norm(a)
    sin = a[2]/norm(a)
    return [cos, sin]
end

function random_start(num_particles::Int64,x_lim::Array{Float64,1},y_lim::Array{Float64,1},vel::Float64)
    size = [abs(x_lim[1]-x_lim[2]),abs(y_lim[1]-y_lim[2])]
    stars = Star[]
    for i in 1:num_particles
        θ = rand()*2*π
        p1 = rand()*2-1
        p2 = rand()*2-1
        star = [Star([p1*size[1]+x_lim[1],p2*size[1]+y_lim[1]],[vel*cos(θ),vel*sin(θ)],10^4)]
        append!(stars,star)
    end
    return stars
end

function cloud_start(num_particles::Int64, center::Array{Float64,1}, c_vel::Array{Float64,1})
    particles = Star[]
    for i = 1:num_particles
        θ = 2π*rand()
        R = randexp()/ℯ * 10.
        z = (rand() > 0.5 ? 1 : -1) * (randexp()/(10ℯ))
        v = √(R*num_particles/1e6)
        mass = 1e4
        push!(particles, Star(center+[R*cos(θ), R*sin(θ)], [v*sin(θ), -v*cos(θ)] + c_vel, mass))
    end
    particles
end

function galaxy_start(num::Int64,center::Array{Float64,1},c_vel::Array{Float64,1},maxR::Float64)
    stars = Star[]
    for i in 1:num
        θ = rand()*2*π
        R = abs(randn()*maxR/2) + maxR/2
        pos = center + [R*cos(θ),R*sin(θ)]
        vel = [sin(θ),-cos(θ)]*29791.032 + c_vel
        append!(stars,[Star(pos,vel,6*10^4)])
    end
    append!(stars,[Star(center,c_vel,1.989*10^10)])
    stars
end

function mean(a::Number,b::Number)
    (b+a)/2
end

function build_animation(history::Array{Float64,3},x_lim::Union{Nothing,Array{Float64,1}},y_lim::Union{Nothing,Array{Float64,1}};df::Int64=1,clash::Bool=false,N::Int64=0,time_unit::String="h")
    println("Building Animation...")
    T = size(history)[1]
    
    anim = @animate for t = tqdm(1:df:T)
        title = time_unit == "h" ? "Time elapsed: $(t÷(24*30*12))y $(t÷(24*30)%12)m" : "Time elapsed: $(t÷(30*12))y $(t÷(30)%12)m"
        if ! clash
            scatter(history[t,1:end-1,1], history[t,1:end-1,2], title=title, legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markersize=3)
            scatter!(history[t,end:end,1], history[t,end:end,2], legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markersize=5)
        else
            scatter(history[t,1:N,1], history[t,1:N,2], title=title, legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markercolor=:deepskyblue3, markersize=3)
            scatter!(history[t,N+1:N+1,1],history[t,N+1:N+1,2], legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markercolor=:deepskyblue3, markersize=5)
            scatter!(history[t,N+2:end-1,1], history[t,N+2:end-1,2], legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markercolor=:orange, markersize=3)
            scatter!(history[t,end:end,1], history[t,end:end,2], legend=false, size=[500,500], xlim=x_lim, ylim=y_lim, markercolor=:orange, markersize=5)
        end
    end
    println("Done!")
    return anim
end

function moveStar(force::Array{Float64,1},star::Star,time::Float64,spaceScale::Int64;sun::Bool=false)
    if sun
        return star
    else
        a = force / star.m * 10^-11 # m/s^2
        v_x = a[1] * time + star.v[1]
        v_y = a[2] * time + star.v[2]
        s_x = a[1]/2 * time^2 + star.v[1] * time + star.s[1] * spaceScale # m
        s_y = a[2]/2 * time^2 + star.v[2] * time + star.s[2] * spaceScale # m 
        return Star([s_x,s_y]/spaceScale,[v_x,v_y],star.m)
    end
end

function collapse(star1::Star,star2::Star)
    if star2.m >= star1.m
        S = star1.s
        V = star1.v
    end
    S = star1.m >= star2.m ? star1.s : star2.s
    V = star1.m >= star2.m ? star1.v : star2.v
    M = star1.m + star2.m
    Star(S,V,M)
end

function mergeStarArrays(array1::Array{Star,1},array2::Array{Star,1},addEarthSun::Bool)
    array = Star[]
    append!(array,array1)
    append!(array,array2)
    array
end

function addEarthSun(array::Array{Star,1},center::Array{Float64,1})
    append!(array,[Star([center[1],center[2]+14.96],[29791.032,0.],5.972*10^4)])
    append!(array,[Star(center,[0.,0.],1.989*10^10)])
    array
end


function viewStart(stars::Array{Star,1})
    position = zeros(length(stars),2)
    for i in 1:length(stars)
        position[i,:] = stars[i].s
    end
    scatter(position[:,1],position[:,2],size=[500,500],legend=false,markersize=5)
end

function viewStart!(stars::Array{Star,1})
    position = zeros(length(stars),2)
    for i in 1:length(stars)
        position[i,:] = stars[i].s
    end
    scatter!(position[:,1],position[:,2],size=[500,500],legend=false)
end