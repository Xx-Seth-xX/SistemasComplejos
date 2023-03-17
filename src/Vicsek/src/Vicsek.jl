module Vicsek
using SCUtils
using Parameters
using StaticArrays
using ProgressMeter

struct Bird
    position::SVector{2, Float64}
    celerity::Float64
    angle::Float64
    velocity::SVector{2, Float64}
    Bird(pos, celerity, angle, velocity) = new(SVector{2}(pos), celerity, angle, velocity)
    Bird(pos, celerity, angle) = new(SVector{2}(pos), celerity, angle, SVector{2}(cos(angle), sin(angle)) * celerity)
    Bird(pos::Vector{Float64}, celerity, angle) = new(SVector{2}(pos), celerity, angle, SVector{2}(cos(angle), sin(angle)) * celerity)
end
angle(bird::Bird) = bird.angle
angle() = zero(Float64)
position(bird::Bird) = bird.position
position() = zero(SVector{2, Float64})
velocity(bird::Bird) = bird.velocity
velocity() = zero(SVector{2, Float64})

@with_kw struct SimulationParameters
    r::Float64
    L::Float64
    N::Float64
    η::Float64
    celerity::Float64
    Δt::Int = 0
    duration::Int
end

generate_random_flock(sim::SimulationParameters) = generate_random_flock(sim.N, sim.L, sim.celerity)
function generate_random_flock(N, L, celerity)
    return [Bird(rand(2) .* L, celerity, rand()*2*π) for _ = 1:N]
end

function calc_new_angle(near_birds::Vector{Bird}, η::Real)
    mean_θ = SCUtils.average(angle, near_birds)
    return mean_θ + (rand()-0.5)* η
end
function calc_new_position(bird::Bird, L::Real)
    return mod.((bird.position + bird.velocity), L)
end

function new_birds(birds::Vector{Bird}, nn_list::Vector{Vector{Int}}, η::Real, L::Real)
    new_birds = Vector{Bird}(undef, length(birds))
    for i in eachindex(nn_list)
        this_p_nn_list = nn_list[i]
        new_angle = calc_new_angle(birds[this_p_nn_list], η)
        new_position = calc_new_position(birds[i], L) 
        new_celerity = birds[i].celerity
        new_velocity = SVector{2}(cos(new_angle), sin(new_angle)) .* new_celerity
        new_birds[i] = Bird(new_position, new_celerity, new_angle, new_velocity)
    end
    return new_birds
end

function step_sim!(flock::Vector{Bird}, sim::SimulationParameters)
    #First we generate the nearest neighbours list
    nn_list = SCUtils.find_nn(position, flock, sim.r, sim.L)
    flock .= new_birds(flock, nn_list, sim.η, sim.L)
end

function calculate_velocity_parameter(flock::Vector{Bird}, sim::SimulationParameters)
    return (sum(velocity, flock) |> SCUtils.norm2 |> sqrt) / (sim.N * sim.celerity)
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
