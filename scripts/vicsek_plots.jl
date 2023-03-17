using DrWatson
@quickactivate
using Vicsek
using GLMakie

import Vicsek.SimulationParameters
DrWatson.allaccess(::SimulationParameters) = (:L, :r, :N, :η, :celerity)
function execute_sim_dict(sim::SimulationParameters) 
    time, va = Vicsek.execute_sim(sim)
    return @strdict(time, va)
end
function execute_sim_keep_frames_dict(sim::SimulationParameters)
    time, flock = Vicsek.execute_sim_keep_frames(sim)
    return @strdict(time, flock)
end

function plot_va_vs_time_at_various_η(;N, L, r, celerity, η)
    duration = 1000
    Δt = 1
    folder = datadir("preliminary_analysis")
    configs = merge(@dict(L, r,  N, celerity, duration, Δt), @dict(η)) |> dict_list
    fig = Figure()
    ax = Axis(fig[1,1])
    ylims!(ax, 0,1)
    for config in configs
        data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
        lines!(ax, data["time"], data["va"], label = "η = $(config[:η])")
    end
    axislegend(ax)
    safesave(plotsdir("preliminary_analysis", savename(@dict(L, r,  N, celerity), "png")),
               fig)
    return fig
end

function animate_simulation(; N, L, r, celerity, η)
    framerate = 60
    duration = 5000
    config = @dict N L r celerity η duration
    fig = Figure()
    ax = Axis(fig[1,1])
    counter = Observable(1)
    folder = datadir("preliminary_analysis", "flock")

    data = produce_or_load(execute_sim_keep_frames_dict, SimulationParameters(;config...), folder)[1]
    positions = @lift(data["flock"][$counter] .|> Vicsek.position)
    directions = @lift(data["flock"][$counter] .|> Vicsek.velocity  )
    time = @lift(data["time"][$counter])

    xs = @lift(getindex.($positions,1))
    ys = @lift(getindex.($positions,2))
    us = @lift(getindex.($directions,1))
    vs = @lift(getindex.($directions,2))

    arrows!(ax, xs, ys, us, vs, color = :black, size = 1)

    record(fig, plotsdir("video", savename(SimulationParameters(;config...), "mp4")), 1:(duration + 1), framerate = framerate) do i
        counter[] = i
    end
end
