module AntColony

export aco

using Statistics


"""
    aco(dist_mat; start_node=nothing, end_node=nothing; is_tour=false, beta=2, rho=0.1,
        q=0.1, Q=1.0, tau_min=1.0, tau_max=10.0, max_iter=20, reset_iter=10,
        top_perc_ants=0.05, verbose=false)

Ant Colony Optimization for directed graphs based on Max-Min Ant System, Elitist Ant System and Ant Colony System.

# Arguments
- `dist_mat`: must be a square distance matrix (NxN) and edges are accessed like this `dist_mat[to, from]`
- `start_node`: must be an integer <= N or `nothing` (both nodes or none must be specified)
- `end_node`: must be an integer <= N or `nothing` (both nodes or none must be specified)
- `is_tour`: whether path should start and end with same node (will set `end_node = start_node`)
- `beta`: influences how important distance is for the ants decision making
- `rho`: pheromones evaporation percentage for each iteration (`0.0 <= rho <= 1.0`)
- `q`: percentage of roulette wheel style decisions of ants
- `Q`: pheromone deposit factor
- `tau_min`: lower bound for pheromone levels for any edge
- `tau_max`: upper bound for pheremone levels for any edge (initial level of pheromones)
- `max_iter`: number of iterations run for the algorithm
- `reset_iter`: number of iterations after which the pheromone levels are reset to `tau_max`
- `top_perc_ants`: the percentage of top perfoming ants which deposit pheromones on their trails
- `verbose`: print to console whenever a better path was found

# Examples
## Find a tour
Find any path with the same start and end node i.e. a tour. Note that the start node only appears once!

```julia-repl
julia> distances = rand(5, 5);
julia> aco(distances, is_tour = true)
5-element Array{Int64,1}:
 4
 3
 5
 1
 2
````

## Find a path
Finda a specific path from node 1 to node 5

```julia-repl
julia> distances = rand(5, 5);
julia> aco(distances, start_node = 1, end_node = 5)
5-element Array{Int64,1}:
 1
 2
 4
 3
 5
````
"""
function aco(
        dist_mat::AbstractMatrix{<:Number};
        start_node::Union{Nothing, Int} = nothing,
        end_node::Union{Nothing, Int} = nothing,
        is_tour = false,
        beta = 2,
        rho = 0.1,
        q = 0.1,
        Q = 1.0,
        tau_min = 1.0,
        tau_max = 10.0,
        max_iter = 20,
        reset_iter = 10,
        top_perc_ants = 0.05,
        verbose = false
    )

    @assert isnothing(start_node) && isnothing(end_node) || !isnothing(start_node) && !isnothing(end_node)
    @assert if is_tour && (!isnothing(start_node) || !isnothing(end_node)); start_node != end_node else true end

    _, n_nodes = size(dist_mat)
    ants = trunc(Int, n_nodes)
    top_k_ants = max(1, trunc(Int, top_perc_ants * ants))

    # how many factors an edge is better than the expected edge (median) to that node
    η = (median(dist_mat, dims=1)' ./ dist_mat) .^ beta
    τ = fill(tau_max, n_nodes, n_nodes)

    no_improv = 0
    best_path = nothing
    best_cost = Inf

    for i in 1:max_iter
        # probability matrix for the ants decision making
        P = τ .* η
        solutions = []
        _start_node, _end_node = start_node, end_node

        for _ in 1:ants
            # select random start/end nodes
            if isnothing(start_node)
                _start_node = rand(1:n_nodes)
                if is_tour
                    _end_node = _start_node
                end
            end
            path = travel(P, q, _start_node, _end_node)
            cost = sum(edge_distances(dist_mat, path, is_tour))
            push!(solutions, (cost, path))
        end

        # update if better solution found
        sort!(solutions)
        best_local_cost, best_local_path = solutions[1]
        if best_local_cost < best_cost
            if verbose
                println("Better solution found with cost $(best_local_cost) at iteration $(i)")
            end
            best_path = best_local_path
            best_cost = best_local_cost
            no_improv = 0
        else
            no_improv += 1
        end

        # deposit pheromones
        for (cost, path) in solutions[1:top_k_ants]
            Δτ = Q / cost
            for (from, to) in edges(path, is_tour)
                τ[to, from] += Δτ
            end
        end

        # reset when stagnation
        if no_improv > reset_iter
            no_improv = 0
            τ = fill(tau_max, n_nodes, n_nodes)
        else
            # evaporate pheromones
            τ = max.(tau_min, min.(tau_max, τ * (1 - rho)))
        end

    end

    best_path
end


"""
    travel(P, q, start_node, end_node)

Travel through graph based on probability `P`
"""
function travel(
        P::Matrix{<:Number},
        q::Float64,
        start_node::Int,
        end_node::Union{Nothing, Int} = nothing
    )

    not_visited = trues(size(P)[1])
    n_nodes = count(not_visited)
    not_visited[start_node] = false
    if !isnothing(end_node)
        not_visited[end_node] = false
    end
    path = Array{Int}(undef, n_nodes)
    path[1] = start_node

    node = start_node

    for i in 2:n_nodes
        if i == n_nodes && !isnothing(end_node) && start_node != end_node
            not_visited[end_node] = true
        end

        neighbor_idx = if rand() < q
            sample(1:count(not_visited), P[not_visited, node])
        else
            argmax(P[not_visited, node])
        end
        next_node = findall(not_visited)[neighbor_idx]

        not_visited[next_node] = false
        path[i] = next_node
        node = next_node
    end
    path
end


"""
    edges(path, is_tour)

Generate node pairs representing edges in `path`
"""
function edges(path::AbstractArray{Int}, is_tour::Bool)
    path_len = length(path)
    n_edges = if is_tour path_len else path_len - 1 end
    (
        (path[i], path[mod1(i + 1, path_len)])
        for i in 1:n_edges
    )
end


"""
    edge_distances(dist_mat, path, is_tour)

Generate a distances for each edge in `path`
"""
function edge_distances(
        dist_mat::AbstractMatrix{<:Number},
        path::AbstractArray{Int},
        is_tour::Bool
    )
    (dist_mat[to, from] for (from, to) in edges(path, is_tour))
end


"""
    sample(x, weights)

Sample from `x` with a custom distribution defined by `weights`
"""
function sample(x, weights)
    ecdf = cumsum(weights)
    p = rand() * ecdf[end]
    for i in 1:length(ecdf)
        if p <= ecdf[i]
            return x[i]
        end
    end
end


end
