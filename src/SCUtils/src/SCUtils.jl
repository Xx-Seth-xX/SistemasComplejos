module SCUtils
using StaticArrays

#Utilities for finding near neighbours
struct CellList
    cells::Dict{SVector{2, Int}, Vector{Int}}
end
function CellList(p::Vector{T}, r::Real) where T
    return CellList(identity, p, r)
end
function CellList(foo::Function, p::Vector{T}, r::Real) where T
    #We make sure r is > 0
    @assert r > 0
    data = Dict{SVector{2, Int}, Vector{Int}}()
    for i in eachindex(p)
        #We calculate the associated cell for that particle
        cell = Int.(div.(foo(p[i]),r,RoundUp))
        #If the cell is already registered in the data we add the index of that particle to it's vector, otherwise we push the cell to data with said index
        if haskey(data, cell)
            push!(data[cell], i)
        else
            data[cell] = [i]
        end
    end

    return CellList(data)
end


find_nn(positions::Vector{T}, r::Real, L::Real) where T = 
    find_nn(identity, positions, r, L)
function find_nn(foo::Function, positions::Vector{T}, r::Real, L::Real) where T 
    find_nn!(Vector{Vector{Int}}(undef, length(positions)), foo, positions, r, L)
end
find_nn!(nn_list::Vector{Vector{Int}}, positions::Vector{T}, r::Real, L::Real) where T  = 
    find_nn!(identity, nn_list, positions, CellList(positions, r), float(r), float(L))
function find_nn!(nn_list::Vector{Vector{Int}}, foo::Function, positions::Vector{T}, r::Real, L::Real) where T  
    find_nn!(nn_list, foo,positions, CellList(foo, positions, r), float(r), float(L))
end
function find_nn!(nn_list::Vector{Vector{Int}}, foo::Function,positions::Vector{T}, data::CellList, r::Float64,L::Float64) where T 
    #We first generate the offsets to all adjacent cells
    cells = data.cells
    offsets = [SVector(i,j) for i = -1:1, j=-1:1]

    max_cell_index = Int(div(L, r, RoundUp))

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
                            push!(particles_nn_to_p, alter_p)
                            continue
                        end
                        if norm2((foo(positions[alter_p]) + pos_offset) - foo(positions[p])) < r2
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
            if alter_i == 
                push!(particles_nn_to_i, alter_i)
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


function average(values::Vector{T}) where T <: Real
    N = length(values)
    return sum(values, init = zero(T)) / N
end
function average(f::Function, values::Vector{T}) where T <: Real
    N = length(values)
    return sum(values, init = zero(T)) / N
end
function average(f::Function, values::Vector{T}) where T 
    N = length(values)
    return sum((x) -> f(x), values, init = f()) /N
end

end # module SCUtils
