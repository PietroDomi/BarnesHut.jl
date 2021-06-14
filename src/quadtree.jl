using Parameters

@with_kw mutable struct Node2D
    parent = nothing
    children = Dict{String,Union{Nothing,Node2D}}([("nw",nothing),
                                                   ("ne",nothing),
                                                   ("sw",nothing),
                                                   ("se",nothing)])
    stars::Array{Star,1} = []
    regCenter::Array{Float64,1} = [0.,0.]
    x_lim::Array{Float64,1} = [0.,0.]
    y_lim::Array{Float64,1} = [0.,0.]
    n::Int64 = 0
    centerOfMass::Array{Float64,1} = [0.,0.]
    totalMass::Float64 = 0.
end

function getQuadrant(center::Array{Float64,1},point::Array{Float64,1})
    if point[1] <= center[1]
        if point[2] <= center[2]
            return "sw"
        else
            return "nw"
        end
    else
        if point[2] > center[2]
            return "ne"
        else
            return "se"
        end
    end
end

function createChild(node::Union{Nothing,Node2D},parent::Node2D,star::Star,
                     center::Array{Float64,1},x_lim::Array{Float64,1},y_lim::Array{Float64,1},level::Int64)
    if isnothing(node)
        node = Node2D(parent=parent,level=level,
                              regCenter=center,x_lim=x_lim,y_lim=y_lim)
        node = updateNode!(node,star)
        return true, node
    else
        node = updateNode!(node,star)
        level += 1
    end
end

function updateNode!(node::Node2D, star::Star)
    node.centerOfMass = [node.totalMass*node.centerOfMass[1] + star.m*star.s[1],
                         node.totalMass*node.centerOfMass[2] + star.m*star.s[2]]
    node.totalMass += star.m
    node.centerOfMass /= node.totalMass
    node.n += 1
    append!(node.stars,[star])
    return node
end

function createRoot(stars::Array{Star,1},x_lim::Array{Float64,1},y_lim::Array{Float64,1})
    center = [mean(x_lim[1],x_lim[2]), mean(y_lim[1],y_lim[2])]
    node = Node2D(regCenter=center,n=length(stars),x_lim=x_lim,y_lim=y_lim)
    node.stars = copy(stars)
    cm = [0.,0.]
    tm = 0.
    for s in node.stars
        cm = [cm[1] + s.s[1]*s.m, cm[2] + s.s[2]*s.m]
        tm += s.m
    end
    cm /= tm
    node.centerOfMass = cm
    node.totalMass = tm
    return node
end

function getCenterLimits(node::Node2D,Q::String)
    center = [mean(node.x_lim[1],node.regCenter[1]),mean(node.y_lim[1],node.regCenter[2])]
    x_lim = [node.x_lim[1],node.regCenter[1]]
    y_lim = [node.y_lim[1],node.regCenter[2]]
    if Q == "nw"
        center[2] = mean(node.regCenter[2],node.y_lim[2])
        y_lim = [node.regCenter[2],node.y_lim[2]]
    elseif Q == "ne"
        center = [mean(node.regCenter[1],node.x_lim[2]),mean(node.regCenter[2],node.y_lim[2])]
        x_lim = [node.regCenter[1],node.x_lim[2]]
        y_lim = [node.regCenter[2],node.y_lim[2]]
    elseif Q == "se"
        center[1] = mean(node.regCenter[2],node.x_lim[2])
        x_lim = [node.regCenter[2],node.x_lim[2]]
    end
    return center, x_lim, y_lim
end

function buildQTree(root::Union{Nothing,Node2D},stars::Union{Nothing,Array{Star,1}},x_lim::Array{Float64,1},y_lim::Array{Float64,1})
    # regCenter = [mean(x_lim[1],x_lim[2]), mean(y_lim[1],y_lim[2])]
    if isnothing(root)
        root = createRoot(stars,x_lim,y_lim)        
    end
#     println(root,stars)
    for s in stars
#         println(i)
        Q = getQuadrant(root.regCenter,s.s)
        if isnothing(root.children[Q])
#             println("bottom")
            center, x_lim, y_lim = getCenterLimits(root,Q)
            root.children[Q] = Node2D(parent=root, regCenter=center,x_lim=x_lim,y_lim=y_lim)
            root.children[Q] = updateNode!(root.children[Q],s)
        else
#             println("going")
            root.children[Q] = updateNode!(root.children[Q],s)
            if root.children[Q].n == 2
                stars_ = root.children[Q].stars
            else
                stars_ = [s]
            end
            root.children[Q] = buildQTree(root.children[Q], stars_, root.children[Q].x_lim, root.children[Q].y_lim)
        end
    end
    return root
end

function computeForceTree(node::Node2D,star::Star,theta::Float64)
    F_s = [0.,0.]
    s = abs(node.x_lim[2]-node.x_lim[1])
    d = distance(star.s,node.centerOfMass)
    if node.n == 1
        if star ∉ node.stars
            f = newton(star,node.stars[1])
            dir = node.stars[1]-star
            cos0, sin0 = get_cos_sin(dir)
            F_s += [f * cos0, f * sin0] # * 10^9
        else
            F_s += [0.,0.]
        end
    elseif s/d < theta
        cm = Star(node.centerOfMass,[0.,0.],node.totalMass)
        f = newton(star,cm)
        dir = cm-star
        cos0, sin0 = get_cos_sin(dir)
        F_s += [f * cos0, f * sin0] # * 10^9
    else
        for k in keys(node.children)
            if ! isnothing(node.children[k])
                F_s += computeForceTree(node.children[k],star,theta)
            end
        end
    end
    return F_s
end


function forcesVector(root::Node2D, theta::Float64)
    F = zeros(length(root.stars),2)
    # This process could be parallelized
    for i in 1:length(root.stars)
        F[i,:] = computeForceTree(root,root.stars[i],theta)
    end
    return F
end


function oneStepTree(stars::Array{Star,1},delta::Float64,theta::Float64,spaceScale::Int64,eps::Float64)
#     println("Finding borders")
    x_lim = [Inf,-Inf]
    y_lim = [Inf,-Inf]
    for s in stars
        x_lim = [min(x_lim[1],s.s[1]),max(x_lim[2],s.s[1])]
        y_lim = [min(y_lim[1],s.s[2]),max(y_lim[2],s.s[2])]
    end
    x_lim = [x_lim[1]-eps,x_lim[2]+eps]
    y_lim = [y_lim[1]-eps,y_lim[2]+eps]
#     println("Building Tree")
    root = buildQTree(nothing,stars,x_lim,y_lim)
#     println("Computing Forces")
    F = forcesVector(root,theta)
    new_stars = copy(stars)
#     println("Moving stars")
    for i in 1:length(stars)
        #TODO: this could be improved exploiting matrix multiplications
         new_stars[i] = moveStar(F[i,:],stars[i],delta,spaceScale)
    end
    return new_stars
end

function simulationTree(stars::Array{Star,1},time::Int64,delta::Float64,theta::Float64,plotStart::Bool)
    position = zeros(length(stars),2)
    for i in 1:length(stars)
        position[i,:] = stars[i].s
    end
    if plotStart
        scatter(position[:,1],position[:,2])
    end
    history = [position]
    stars_ = copy(stars)
    println("Beginning BH simulation...")
    for t in 1:time 
        stars_ = oneStepTree(stars_,delta,theta,10^10,1.)
        if t % (time÷20) == 1
            println("$t out of $time")
        end
        position = zeros(length(stars_),2)
        for i in 1:length(stars_)
            position[i,:] = stars_[i].s
        end
        append!(history, [position])
    end
    return history
end