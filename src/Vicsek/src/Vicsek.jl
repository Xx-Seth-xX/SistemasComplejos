module Vicsek
using SCUtils
using Parameters
using StaticArrays
using ProgressMeter
using Statistics

struct Bird
    position::SVector{2, Float64}
    celerity::Float64
    velocity::SVector{2, Float64}
    Bird(position::SVector{2, Float64}, celerity::Float64, velocity::SVector{2, Float64}) = new(position, celerity, velocity)
    Bird(pos::Vector{Float64}, celerity::Float64, angle::Float64) = new(SVector{2}(pos), celerity, SVector{2}(cos(angle), sin(angle)) .* celerity)
end
position(bird::Bird) = bird.position
position() = zero(SVector{2, Float64})
velocity(bird::Bird) = bird.velocity
velocity() = zero(SVector{2, Float64})

@with_kw struct SimulationParameters
    r::Float64 = 1.0
    L::Float64
    N::Float64
    η::Float64
    ϕ::Float64 = π
    celerity::Float64 = 0.03
    Δt::Int = 0
    duration::Int
end

generate_random_flock(sim::SimulationParameters) = generate_random_flock(sim.N, sim.L, sim.celerity)
function generate_random_flock(N, L, celerity)
    return [Bird(rand(2) .* L, celerity, rand()*2*π) for _ = 1:N]
end

function calc_new_angle(near_birds::Vector{Bird},this::Bird, η::Real, ϕ::Float64, r::Float64, L::Float64)
    birds_in_cone = Iterators.filter((other_bird) -> check_if_bird_in_vision_cone(this, other_bird, ϕ, r, L), Iterators.Stateful(near_birds))
    mean_v = mean((bird) -> bird.velocity, birds_in_cone)
    mean_θ = atan(mean_v.y, mean_v.x)
    return mean_θ + (rand()-0.5)* η
end
function calc_new_position(bird::Bird, L::Real)
    return mod.((bird.position + bird.velocity), L)
end

function check_if_bird_in_vision_cone(this::Bird, other::Bird, ϕ::Float64, r::Float64, L::Float64)
    # Angle between line pointing from this to other and horizontal
    diff_x = other.position[1] - this.position[1]
    diff_y = other.position[2] - this.position[2]
    if diff_y == 0 && diff_x == 0
        return true
    end
    if diff_x < -r 
        diff_x += L
    elseif diff_x > r
        diff_x -= L
    end
    if diff_y < -r 
        diff_y += L
    elseif diff_y > r
        diff_y -= L
    end
    α = atan(diff_y, diff_x)
    # Angle between this bird pointing vector and horizontal
    θ = atan(this.velocity[2], this.velocity[1])
    return abs(α - θ) ≤ ϕ
end

function new_birds(birds::Vector{Bird}, nn_list::Vector{Vector{Int}}, η::Real, L::Real, r::Float64, ϕ::Float64)
    new_birds = Vector{Bird}(undef, length(birds))
    Threads.@threads for i in eachindex(nn_list)
        this_p_nn_list = nn_list[i]
        new_angle = calc_new_angle(birds[this_p_nn_list],birds[i], η, ϕ, r, L)
        new_position = calc_new_position(birds[i], L) 
        new_celerity = birds[i].celerity
        new_velocity = SVector{2}(cos(new_angle), sin(new_angle)) .* new_celerity
        new_birds[i] = Bird(new_position, new_celerity, new_velocity)
    end
    return new_birds
end

function step_sim!(flock::Vector{Bird}, sim::SimulationParameters)
    #First we generate the nearest neighbours list
    nn_list = SCUtils.find_nn(position, flock, sim.r, sim.L)
    flock .= new_birds(flock, nn_list, sim.η, sim.L, sim.r, sim.ϕ)
end

function calculate_velocity_parameter(flock::Vector{Bird}, sim::SimulationParameters)
    return (sum(velocity, flock) |> SCUtils.norm2) / (sim.N * sim.celerity)^2
end

function execute_sim(sim::SimulationParameters)
    flock = generate_random_flock(sim.N, sim.L, sim.celerity)
    
    
    #Saving timeseries
    time = Int[]
    va = Float64[]
    
    push!(time, 0)
    push!(va, calculate_velocity_parameter(flock, sim))
    last_saved_time = 0

    @showprogress 1 for step_number in 1:sim.duration
        step_sim!(flock, sim)
        if(step_number - last_saved_time) ≥ sim.Δt
            last_saved_time = step_number
            push!(time, step_number)
            push!(va, calculate_velocity_parameter(flock, sim))
        end
    end
    return time, va
end

function execute_sim_keep_frames(sim::SimulationParameters)
    flock = generate_random_flock(sim.N, sim.L, sim.celerity)
    
    
    #Saving timeseries
    time = Int[]
    flock_v = Vector{Bird}[]
    
    push!(time, 0)
    push!(flock_v, flock)
    last_saved_time = 0

    @showprogress 1 for step_number in 1:sim.duration
        step_sim!(flock, sim)
        if(step_number - last_saved_time) ≥ sim.Δt
            last_saved_time = step_number
            push!(time, step_number)
            push!(flock_v, copy(flock))
        end
    end
    return time, flock_v
end
end # module Vicsek
