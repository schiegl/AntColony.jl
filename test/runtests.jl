using Test

include("../src/AntColony.jl")
using .AntColony: aco, travel


@testset "Travel a specific path" begin
    n = 30
    P = rand(n, n)
    for start_node in 1:n, end_node in 1:n
        if start_node != end_node
            path = travel(P, 0.1, start_node, end_node)
            @test path[1] == start_node
            @test path[end] == end_node
            @test sort(path) == 1:n
            @test setdiff(Set(path), Set(1:n)) |> isempty
        end
    end
end


@testset "Travel any path" begin
    n = 30
    P = rand(n, n)
    for start_node in 1:n
        path = travel(P, 0.1, start_node, nothing)
        @test path[1] == start_node
        @test sort(path) == 1:n
        @test setdiff(Set(path), Set(1:n)) |> isempty
    end
end


@testset "Travel a tour" begin
    n = 30
    P = rand(n, n)
    for start_node in 1:n
        path = travel(P, 0.1, start_node, start_node)
        @test path[1] == start_node
        @test sort(path) == 1:n
        @test setdiff(Set(path), Set(1:n)) |> isempty
    end
end


@testset "ACO constructs a specific path" begin
    n = 20
    dist_mat = rand(n, n)
    for start_node in 1:n, end_node in 1:n
        if start_node != end_node
            path = aco(dist_mat, start_node=start_node, end_node=end_node, is_tour=false, max_iter=2)
            @test path[1] == start_node
            @test path[end] == end_node
            @test sort(path) == 1:n
            @test setdiff(Set(path), Set(1:n)) |> isempty
        end
    end
end


@testset "ACO constructs any path" begin
    n = 20
    dist_mat = rand(n, n)
    for _ in 1:n
        path = aco(dist_mat, is_tour=false, max_iter=2)
        @test sort(path) == 1:n
        @test setdiff(Set(path), Set(1:n)) |> isempty
    end
end


@testset "ACO constructs a tour" begin
    n = 20
    dist_mat = rand(n, n)
    for start_node in 1:n
        path = aco(dist_mat, is_tour=true, max_iter=2)
        @test sort(path) == 1:n
        @test setdiff(Set(path), Set(1:n)) |> isempty
    end
end
