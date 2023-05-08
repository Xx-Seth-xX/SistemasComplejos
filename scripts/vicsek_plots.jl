using DrWatson
@quickactivate
using Vicsek
using GLMakie
using Statistics
using StatsBase
using CSV, CSVFiles
using DataFrames

import Vicsek.SimulationParameters
function execute_sim_dict(sim::SimulationParameters) 
    time, va = Vicsek.execute_sim(sim)
    return @strdict(time, va)
end
function execute_sim_keep_frames_dict(sim::SimulationParameters)
    time, flock = Vicsek.execute_sim_keep_frames(sim)
    return @strdict(time, flock)
end

function plot_va_vs_time_at_various_η(;N, L, r, celerity, η)
    duration = 10_000
    folder = datadir("preliminary_analysis")
    configs = merge(@dict(L, r,  N, celerity, duration), @dict(η)) |> dict_list
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

function calculate_correlation_time(; N, L)
    η = 3
    duration = 1_000_000
    Δt = 1 #We need to keep every frame so we can calculate correlation time
    folder = datadir("preliminary_analysis", "correlation_time")
    config = @dict(N, L, duration, Δt, η)
    println("Calculating correlation time for N = $(config[:N]), L = $(config[:L]), η = $(config[:η])")
    data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]

    #We get va from data
    va = data["va"]

    #We ignore the first 10_000 va points so we can get rid of the transient behaviour
    va = va[10_000:end]

    #We calculate the autocorrelation function
    acf = autocor(va)

    #We calculate the correlation time as the first data point that goes below 1/e
    τ = findfirst(x -> x < 1/exp(1), acf)
    println("Correlation time for given config is: $(τ)")
    return τ
end

function plot_correlation_time(; N, L, η)
    duration = 1_000_000
    Δt = 1 #We need to keep every frame so we can calculate correlation time
    folder = datadir("preliminary_analysis")
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Tiempo", ylabel = "Autocorrelación")
    τ = Float64
    for (N, L, η) in zip(N, L, η)
        config = @dict(N, L, duration, Δt, η)
        data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]

        #We get va from data
        va = data["va"]

        #We ignore the first 500 va points so we can get rid of the transient behaviour
        va = va[100:end]

        #We calculate the autocorrelation function
        acf = autocor(va, 0:500)

        println("Correlation time for config: $(config) is: $(findfirst(x -> x < 1/exp(1), acf))")
        
        #We plot the autocorrelation function
        lines!(ax, 1:length(acf), acf, label = L"η = %$(config[:η]),\; N = %$(config[:N]),\; L = %$(config[:L])")
    end
    axislegend(ax)
    safesave(plotsdir("preliminary_analysis", "correlation_length.png"), fig)
    return fig
end

function calculate_va_vs_noise(;L, N, tT)
    #We parametrise the noise from 0 to 5 at 0.1 intervals
    L = float(L)
    η = 0.1:0.1:5
    number_of_data_points = 500
    for (L, N, tT) in zip(L, N, tT)
        va_mean = Float64[]
        va_error = Float64[]
        va_var = Float64[]
        mvar_error = Float64[]
        folder = datadir("raw_data")
        #Adaptative time step so we always use 2 times the correlation time
        Δt = 2*calculate_correlation_time(;N = N, L = L)
        for η in η

            #The duration will be the number of data points times the time step plus 100 to eliminate the transient behaviour
            duration = number_of_data_points * Δt + tT
            config = @dict(L, N, duration, Δt, η)
            data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
            va = data["va"]
            @show mean(va[end-number_of_data_points:end])
            @show std(va[end-number_of_data_points:end])*2
            #We calculate the mean of the va for each noise value
            va = va[end-number_of_data_points:end]
            @show mean(va)
            @show std(va)*2 / (sqrt(length(va) - 1))
            push!(va_mean, mean(va))
            push!(va_var, std(va, corrected = false))
            push!(va_error, std(va)*2 / sqrt(length(va) - 1)) 
        end
        #we save the data to csv
        safesave(datadir("va_vs_noise", savename(@dict(L, N), "csv")), DataFrame(@strdict(η, va_mean, va_error, va_var)))
    end
end
