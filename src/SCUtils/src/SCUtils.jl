module SCUtils
using Reexport
@reexport using StaticArrays

struct Bird
    position::SVector{3, Float64}
    celerity::Float64
    angle::Float64
end
struct Flock
    birds::Vector{Bird}
end


end # module SCUtils
