module SCUtils
using Reexport
@reexport using Test
@reexport using StaticArrays
@reexport using SparseArrays

#Utilities for finding near neighbours
struct CellList
    cells::Dict{SVector{2, Int}, Vector{Int}}
end
function CellList(p::Vector{SVector{2,Float64}}, r::Float64, L::Float64)
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

function find_nn(positions::Vector{SVector{2, Float64}}, r::Real, L::Real) 
    find_nn(positions, CellList(positions, r, L), float(r), float(L))
end
function find_nn(positions::Vector{SVector{2, Float64}}, data::CellList, r::Float64,L::Float64)
    #We first generate the offsets to all adjacent cells
    cells = data.cells
    offsets = [SVector(i,j) for i = -1:1, j=-1:1]

    max_cell_index = Int(div(L, r, RoundUp))

    nn_list = [Int[] for _ in positions]

    r2 = r^2
    for (cell, particles) in cells
        for p in particles
            particles_nn_to_p = Int[]
            for offset in offsets
                pos_offset = SVector(0.0,0.0)

                ncell = cell + offset

                if ncell[1] == max_cell_index + 1
                    pos_offset += SVector(L, 0.0)
                    ncell = SVector(1,ncell[2])
                elseif ncell[1] == 0
                    pos_offset += SVector(-L, 0.0)
                    ncell = SVector(max_cell_index,ncell[2])
                end

                if ncell[2] == max_cell_index + 1
                    pos_offset += SVector(0.0, L)
                    ncell = SVector(ncell[1], 1)
                elseif ncell[2] == 0
                    pos_offset += SVector(0.0, -L)
                    ncell = SVector(ncell[1], max_cell_index)
                end

                if haskey(cells, ncell)
                    for alter_p in cells[ncell]
                        if alter_p == p
                            continue
                        end
                        if norm2((positions[alter_p] + pos_offset) - positions[p]) < r2
                            push!(particles_nn_to_p, alter_p)
                        end
                    end
                end
            end
            nn_list[p] = particles_nn_to_p
        end
    end
    return nn_list
end

function norm2(v::SVector{d}) where d
    sum = 0.0
    for i in 1:d
        sum += v[i]^2
    end
    return sum
end

#Just for testing the other function, incredibily innefficient, even for a brute force search function
function find_nn_brute_force(particles::Vector{SVector{2, Float64}}, r, L)
    nn_list = [Int[] for _ in 1:length(particles)]
    r2 = r^2
    aux = (L - r)^2
    offsets = [SVector(i, j) for i = [0,-L,L], j = [0,-L,L]]
    for i in eachindex(particles)
        particles_nn_to_i = Int[]
        for alter_i in eachindex(particles)
            if alter_i == i
                continue
            elseif norm2(particles[i] - particles[alter_i]) < r2              
                push!(particles_nn_to_i, alter_i)
            elseif norm2(particles[i] - particles[alter_i]) > aux
                for offset in offsets
                    if norm2(particles[i] - particles[alter_i] + offset) < r2
                        push!(particles_nn_to_i, alter_i)
                        break
                    end
                end
            end
        end
        nn_list[i] = particles_nn_to_i
    end
    return nn_list
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
    @test CellList(particles,r, 3.0).cells == Dict([SVector(1,1) => [1,2], 
                                                SVector(2,1) => [3],
                                                SVector(2,2) => [4,5],
                                                SVector(3,3) => [6]])

end

end # module SCUtils
