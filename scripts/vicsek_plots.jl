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

function plot_va_vs_time_at_various_η(N, L, r, celerity, η)
    duration = 1000
    Δt = 1
    folder = datadir("preliminary_analysis")
    configs = merge(@dict(L, r,  N, celerity, duration, Δt), @dict(η)) |> dict_list
    fig = Figure()
    ax = fig(ax[1,1])
    for config in configs
        produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)
        ax
    end

end