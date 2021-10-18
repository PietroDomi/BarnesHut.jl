# using Parameters

import Base: show

function empty_children()
    Dict([("nw",nothing),("ne",nothing),("sw",nothing),("se",nothing)])
end

mutable struct Node2D
    parent::Union{Nothing,Node2D}
    children::Dict{String,Union{Nothing,Node2D}}
    stars::Array{Star,1}
    regCenter::Array{Float64,1}
    x_lim::Array{Float64,1}
    y_lim::Array{Float64,1}
    n::Int64
    centerOfMass::Array{Float64,1}
    totalMass::Float64
    Node2D(parent,children,stars,regCenter,x_lim,
            y_lim,n,centerOfMass,totalMass) = 
        new(parent,children,stars,regCenter,x_lim,
            y_lim,n,centerOfMass,totalMass)
end

show(io::IO, n::Node2D) = println(io, "Node with $(n.n) stars, cM=$(n.centerOfMass) and rC=$(n.regCenter)")

Node2D() = Node2D(nothing,empty_children(),[],[0.,0.],[0.,0.],[0.,0.],0,[0.,0.],0.)
Node2D(parent,center::Vector{Float64},x_lim,y_lim) = Node2D(parent,empty_children(),Array{Star}([]),center,x_lim,y_lim,0,[0.,0.],0.)
Node2D(center::Vector{Float64},n,x_lim,y_lim) = Node2D(nothing,empty_children(),[],center,x_lim,y_lim,n,[0.,0.],0.)

function get_quadrant(center::Vector{Float64},point::Vector{Float64})
    return point[1] ≤ center[1] ? (point[2] ≤ center[2] ? "sw" : "nw") :
                                  (point[2] > center[2] ? "ne" : "se")
end


# function createChild(node::Union{Nothing,Node2D},parent::Node2D,star::Star,
#                      center::Array{Float64,1},x_lim::Array{Float64,1},y_lim::Array{Float64,1},level::Int64)
#     if isnothing(node)
#         node = Node2D(parent,center,x_lim,y_lim)
#         node = update_node!(node,star)
#         return true, node
#     else
#         node = update_node!(node,star)
#         level += 1
#     end
# end

function update_node!(node::Node2D, star::Star)
    node.centerOfMass = [node.totalMass*node.centerOfMass[1] + star.m*star.s[1],
                         node.totalMass*node.centerOfMass[2] + star.m*star.s[2]]
    node.totalMass += star.m
    node.centerOfMass /= node.totalMass
    node.n += 1
    append!(node.stars,[star])
    return node
end

function create_root(stars::Array{Star,1},x_lim::Array{Float64,1},y_lim::Array{Float64,1})
    center = [mean(x_lim[1],x_lim[2]), mean(y_lim[1],y_lim[2])]
    node = Node2D(center,length(stars),x_lim,y_lim)
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

function get_center_limits(node::Node2D,Q::String)
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

function build_qtree(root::Union{Nothing,Node2D},stars::Union{Nothing,Array{Star,1}},
                    x_lim::Array{Float64,1},y_lim::Array{Float64,1})
    if isnothing(root)
        root = create_root(stars,x_lim,y_lim)        
    end
    # Loop on every star
    for s in stars
        Q = get_quadrant(root.regCenter,s.s)
        if isnothing(root.children[Q])
            # If the quadrant is still empty, create the node
            center, x_lim, y_lim = get_center_limits(root,Q)
            root.children[Q] = Node2D(root,center,x_lim,y_lim)
            root.children[Q] = update_node!(root.children[Q],s)
        else
            # Otherwise add the star to the node
            root.children[Q] = update_node!(root.children[Q],s)
            if root.children[Q].n == 2
                # If there was already a star, then we need to go down a level and create other nodes
                stars_ = root.children[Q].stars
            else
                # otherwise we just need to sort the new star
                stars_ = [s]
            end
            # Then recursively call the constructor for that node and its corresponding quadrant
            root.children[Q] = build_qtree(root.children[Q], stars_, 
                                          root.children[Q].x_lim, root.children[Q].y_lim)
        end
    end
    # At the end it is enough to return the root, as it contains pointers to all its children
    return root
end

function compute_force_tree(node::Node2D,star::Star,theta::Float64)
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
                F_s += compute_force_tree(node.children[k],star,theta)
            end
        end
    end
    return F_s
end

function forces_vector_update(root::Node2D, theta::Float64,delta::Float64,spaceScale::Int64)
    F = zeros(length(root.stars),2)
    # IDEA: This process could be parallelized
    new_stars = copy(root.stars)
    for i in 1:length(root.stars)
        F[i,:] = compute_force_tree(root,root.stars[i],theta)
        # IDEA: this could be improved exploiting matrix multiplications
        new_stars[i] = move_star(F[i,:],root.stars[i],delta,spaceScale)
    end
    return new_stars
end

function one_step_tree(stars::Array{Star,1},delta::Float64,theta::Float64,spaceScale::Int64,eps::Float64)
    # Finding borders
    x_lim = [Inf,-Inf]
    y_lim = [Inf,-Inf]
    for s in stars
        x_lim = [min(x_lim[1],s.s[1]),max(x_lim[2],s.s[1])]
        y_lim = [min(y_lim[1],s.s[2]),max(y_lim[2],s.s[2])]
    end
    x_lim = [x_lim[1]-eps,x_lim[2]+eps]
    y_lim = [y_lim[1]-eps,y_lim[2]+eps]
    # Building Tree
    root = build_qtree(nothing,stars,x_lim,y_lim)
    # Computing Forces & Moving stars
    new_stars = forces_vector_update(root,theta,delta,spaceScale)
    return new_stars
end

function simulation_tree(stars::Array{Star,1},time::Int64,delta::Float64,theta::Float64;plotStart::Bool=false,timing::Bool=false)
    if plotStart
        scatter(position[:,1],position[:,2])
    end
    stars_ = copy(stars)
    history = zeros(time+1,length(stars_),2)
    for i in 1:length(stars)
        history[1,i,:] = stars[i].s
    end
    timing && (bench_time = Float64[])
    println("Beginning BH simulation...")
    for t in tqdm(1:time) 
        if timing
            res = @timed one_step_tree(stars_,delta,theta,10^10,1.)
            stars_ = res.value
            append!(bench_time,[res.time])
        else
            stars_ = one_step_tree(stars_,delta,theta,10^10,1.)
        end
        for i in 1:length(stars_)
            history[t+1,i,:] = stars_[i].s
        end
    end

    if timing
        return history, bench_time
    else
        return history
    end
end