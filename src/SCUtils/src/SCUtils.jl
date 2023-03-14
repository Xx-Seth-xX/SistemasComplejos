module SCUtils
using Reexport
@reexport using Test
@reexport using StaticArrays
@reexport using SparseArrays

#Utilities for finding near neighbours
struct NNList
    list::Vector{Int64}
    NNList() = new(Int64[])
end
function Base.iterate(nn_list::NNList)
    N = nn_list.list[1]
    return nn_list.list[2:(N + 1)], N + 2
end
function Base.iterate(nn_list::NNList, state::Int64)
    if state > length(nn_list.list)
        return nothing
    end
    N = nn_list.list[state]
    return nn_list.list[(state + 1):(state + N)], state + N + 1
end
function add_list(nn_list::NNList, list::Vector{Int64})
    N = length(list)
    append!(nn_list.list, [N list])
end

struct CellList
    cells::Dict{SVector{2, Int}, Vector{Int}}
end
function CellList(p::Vector{SVector{2,Float64}}, r::Float64)
    #We make sure r is > 0
    @assert r > 0
    #We get the number of particles
    n = length(p)
    data = Dict{SVector{2, Int}, Vector{Int}}()
    for i in eachindex(p)
        #We calculate the associated cell for that particle
        cell = Int.(div.(p[i],r,RoundUp))
        #If the cell is already registered in the data we add the index of that particle to it's vector, otherwise we push the cell to data with said index
        if haskey(data, cell)
            push!(data[cell], i)
        else
            data[cell] = [i]
        end
    end
    return CellList(data)
end

#Tests
@testset "CellList" begin
    particles = SVector{2, Float64}[]
    push!(particles, SVector(0.5,0.5)) # (1,1)
    push!(particles, SVector(0.75,0.75)) #(1,1)
    push!(particles, SVector(1.5,0.5)) #(2,1)
    push!(particles, SVector(1.5,1.5)) #(2,2)
    push!(particles, SVector(1.5,1.6)) #(2,2)
    push!(particles, SVector(2.5,2.7)) #(3,3)
    r = 1.0
    @test CellList(particles,r).cells == Dict([SVector(1,1) => [1,2], 
                                                SVector(2,1) => [3],
                                                SVector(2,2) => [4,5],
                                                SVector(3,3) => [6]])

end

end # module SCUtils
