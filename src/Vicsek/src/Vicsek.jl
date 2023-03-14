module Vicsek
using SCUtils

struct Bird
    position::SVector{3, Float64}
    celerity::Float64
    angle::Float64
end

struct Simulation
    L::Int64
    r::Float64
    N::Float64
    flock::Vector{Bird}
    nn_list::Vector{Vector{Int64}}
    function Simulation(L, r, N)

    end
end



end # module Vicsek
