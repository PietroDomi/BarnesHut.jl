function forceDistanceBrute(stars::Array{Star,1})
    F = zeros(length(stars),length(stars))
    d = zeros(length(stars),length(stars),2)
    k = 1
    for i in 1:length(stars)
        for j in 1:k
            if stars[i]!=stars[j]
                f_i = newton(stars[i],stars[j])
                F[i,j] = f_i
                F[j,i] = f_i
                d_i = stars[j]-stars[i]
                d[i,j,1],d[i,j,2] = d_i
                d[j,i,1],d[j,i,2] = -d_i
            end
        end
        k += 1
    end
    return F, d
end

function onestepBrute(time::Float64,stars::Array{Star,1},spaceScale::Int64)
    new_stars = copy(stars)
    for i in 1:length(stars)
        # compute the net force
        f_i = [0,0]
        for j in 1:length(stars)
            if i != j
                F = newton(stars[j],stars[i])
                d_j = stars[j]-stars[i]
                cos_0j, sin_0j = get_cos_sin(d_j)
                f_j = [F * cos_0j, F * sin_0j]
                f_i += f_j # * 10^9
            end
        end
        # move the body according to the net force
        new_stars[i] = moveStar(f_i,stars[i],time,spaceScale)
    end
    return new_stars
end

function simulationBrute(stars::Array{Star,1},time::Int64,delta::Float64;timing::Bool=false)
    position = zeros(length(stars),2)
    for i in 1:length(stars)
        position[i,:] = [stars[i].s[1],stars[i].s[2]]
    end
    history = [position]
    stars_ = copy(stars)
    timing && (bench_time = Float64[])
    println("Beginning brute force simulation...")
    for t in tqdm(1:time) # one month's minutes
        if timing
            res = @timed onestepBrute(delta,stars_,10^10)
            stars_ = res.value
            append!(bench_time,[res.time])
        else
            stars_ = onestepBrute(delta,stars_,10^10)
        end
        # if t % (time÷10) == 1
        #     # println("$t out of $time")
        # end
        position = zeros(length(stars_),2)
        for i in 1:length(stars_)
            position[i,:] = [stars_[i].s[1],stars_[i].s[2]]
        end
        append!(history, [position])
    end
    if timing
        return history, bench_time
    else
        return history
    end
end