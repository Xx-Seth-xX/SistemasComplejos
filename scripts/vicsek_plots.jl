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
    safesave(plotsdir("preliminary_analysis", savename(@dict(L, r,  N, celerity, duration, Δt), "png")),
               fig)
    return fig
end

