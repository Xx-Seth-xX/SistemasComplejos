using DrWatson
@quickactivate
using Vicsek
using ProgressMeter
using CairoMakie
using Statistics
using StatsBase
using CSV, CSVFiles
using DataFrames
using Colors

update_theme!(fontsize = 30)

import Vicsek.SimulationParameters
function execute_sim_dict(sim::SimulationParameters) 
    time, va = Vicsek.execute_sim(sim)
    return @strdict(time, va)
end
function execute_sim_keep_frames_dict(sim::SimulationParameters)
    time, flock = Vicsek.execute_sim_keep_frames(sim)
    return @strdict(time, flock)
end

function plot_va_vs_time(;N, L, η, ϕ)
    duration = 10_000
    folder = datadir("preliminary_analysis")
    configs = merge(@dict(L, N, duration), @dict(η), @dict(ϕ)) |> dict_list
    fig = Figure()
    ax = Axis(fig[1,1])
    ylims!(ax, 0,1)
    for config in configs
        data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
        lines!(ax, data["time"], data["va"], label = L"η = %$(config[:η]), ϕ = %$(config[:ϕ])")
    end
    axislegend(ax)
    safesave(plotsdir("preliminary_analysis", savename(@dict(L, N), "png")),
               fig)
    return fig
end

function animate_simulation(; N, L, η, ϕ, duration = 20)
    framerate = 30
    duration = framerate * duration
    config = @dict N L duration η ϕ
    fig = Figure()
    ax = Axis(fig[1,1])
    xlims!(ax, 0, L)
    ylims!(ax, 0, L)
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
    function get_n_before_i(v::Vector, n, i)  
        if n > i
            n = i
        end
        vout = reduce(hcat, v[(i - n + 1):i])
        return vout
    end

    #rows are the particles and the columns are the times
    tails = @lift(get_n_before_i(data["flock"], 20, $counter) .|> Vicsek.position)
    
    arrows!(ax, xs, ys, us, vs, color = :black, size = L/2)

    for i in 1:size(tails[])[1]
        tail = @lift(($tails)[i,:])
        @show tail
        scatter!(ax, tail, markersize = 6, markerspace = :pixel, color = RGBA(0,1,0,0.3))
    end

    p = Progress(duration + 1, desc = "Rendering simulation...", dt=1.0, showspeed = true)

    record(fig, plotsdir("video", savename(SimulationParameters(;config...), "mp4")), 1:(duration + 1), framerate = framerate) do i
        update!(p, i)
        counter[] = i
    end
end

function calculate_correlation_time(; N, L)
    η = 2
    duration = 1_000_000
    Δt = 1 #We need to keep every frame so we can calculate correlation time
    folder = datadir("preliminary_analysis", "correlation_time")
    config = @dict(N, L, duration, Δt, η)
    println("Calculating correlation time for N = $(config[:N]), L = $(config[:L]), η = $(config[:η])")
    data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]

    #We get va from data
    va = data["va"]

    #We ignore the first 10_000 va points so we can get rid of the transient behaviour
    va = va[300:end]

    #We calculate the autocorrelation function
    acf = autocor(va, 0:500)
    #We calculate the correlation time as the first data point that goes below 1/e
    τ = findfirst(x -> x < exp(-1), acf)
    println("Correlation time for given config is: $(τ)")
    return τ
end
function calculate_correlation_time_angle(; N, L, ϕ = π, η = 2)
    η = 2
    duration = 1_000_000
    Δt = 1 #We need to keep every frame so we can calculate correlation time
    folder = datadir("preliminary_analysis", "correlation_time")
    config = @dict(N, L, duration, Δt, η, ϕ)
    println("Calculating correlation time for N = $(config[:N]), L = $(config[:L]), η = $(config[:η]), ϕ = $(config[:ϕ])")
    data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]

    #We get va from data
    va = data["va"]

    #We ignore the first 10_000 va points so we can get rid of the transient behaviour
    va = va[300:end]

    #We calculate the autocorrelation function
    acf = autocor(va, 0:500)
    #We calculate the correlation time as the first data point that goes below 1/e
    τ = findfirst(x -> x < exp(-1), acf)
    println("Correlation time for given config is: $(τ)")
    return τ
end

function plot_correlation_time(; N, L, η)
    duration = 1_000_000
    Δt = 1 #We need to keep every frame so we can calculate correlation time
    folder = datadir("preliminary_analysis", "correlation_time")
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Tiempo", ylabel = "Autocorrelación")
    τ = Float64
    for (N, L, η) in zip(N, L, η)
        config = @dict(N, L, duration, Δt, η)
        data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]

        #We get va from data
        va = data["va"]

        #We ignore the first 500 va points so we can get rid of the transient behaviour
        va = va[300:end]

        #We calculate the autocorrelation function
        acf = autocor(va, 0:500)

        corr_time = findfirst(x -> x < exp(-1), acf)

        println("Correlation time for config: $(config) is: $(corr_time)")
        
        #We plot the autocorrelation function
        lines!(ax, 1:length(acf), acf, label = L"η = %$(config[:η]),\; N = %$(config[:N]),\; L = %$(config[:L])")
    end
    axislegend(ax)
    safesave(plotsdir("preliminary_analysis", "correlation_time.png"), fig)
    return fig
end

function calculate_va_vs_noise(;L, N)
    #We parametrise the noise from 0 to 5 at 0.1 intervals
    L = float(L)
    η = 0.2:0.2:5
    number_of_data_points = 3000
    for (L, N) in zip(L, N)
        va_mean = Float64[]
        va_error = Float64[]
        va_var = Float64[]
        folder = datadir("raw_data")
        #Adaptative time step so we always use 2 times the correlation time
        Δt = 2*calculate_correlation_time(;N = N, L = L)
        for η in η

            #The duration will be the number of data points times the time step plus 300 to eliminate the transient behaviour
            duration = number_of_data_points * Δt + 300
            config = @dict(L, N, duration, Δt, η)
            data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
            va = data["va"]
            @show mean(va[end-number_of_data_points:end])
            @show std(va[end-number_of_data_points:end])*2
            #We calculate the mean of the va for each noise value
            va = va[end-number_of_data_points:end]
            @show mean(va)
            @show std(va)*2 / (sqrt(length(va) - 1))
            @show std(va)*2 
            push!(va_mean, mean(va))
            push!(va_var, std(va, corrected = false))
            push!(va_error, std(va)*2 / sqrt(length(va) - 1)) 
        end
        #we save the data to csv
        safesave(datadir("va_vs_noise", savename(@dict(L, N), "csv")), DataFrame(@strdict(η, va_mean, va_error, va_var)))
    end
end
function calculate_va_vs_angle(;L, N, η)
    #We parametrise the noise from 0 to 5 at 0.1 intervals
    L = float(L)
    ϕ = range(0, π, length = 15)
    number_of_data_points = 3000
    for (L, N) in zip(L, N)
        va_mean = Float64[]
        va_error = Float64[]
        va_var = Float64[]
        folder = datadir("raw_data")
        #Adaptative time step so we always use 2 times the correlation time
        Δt = 2*calculate_correlation_time_angle(;N = N, L = L, ϕ = π/2, η = 2)
        for ϕ in ϕ
            #The duration will be the number of data points times the time step plus 300 to eliminate the transient behaviour
            duration = number_of_data_points * Δt + 500
            config = @dict(L, N, duration, Δt, η, ϕ)
            data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
            va = data["va"]
            @show mean(va[end-number_of_data_points:end])
            @show std(va[end-number_of_data_points:end])*2
            #We calculate the mean of the va for each noise value
            va = va[end-number_of_data_points:end]
            @show mean(va)
            @show std(va)*2 / (sqrt(length(va) - 1))
            @show std(va)*2 
            @show ϕ
            push!(va_mean, mean(va))
            push!(va_var, std(va, corrected = false))
            push!(va_error, std(va)*2 / sqrt(length(va) - 1)) 
        end
        #we save the data to csv
        safesave(datadir("va_vs_noise", savename(@dict(L, N, η), "csv")), DataFrame(@strdict(ϕ, va_mean, va_error, va_var)))
    end
end
function calculate_va_vs_density_fixed_size(;L, η)
    L = float(L)
    η = float(η)
    number_of_data_points = 3000
    for (L, η) in zip(L, η)
        va_mean = Float64[]
        va_error = Float64[]
        va_var = Float64[]
        folder = datadir("raw_data")
        #Adaptative time step so we always use 2 times the correlation time
        Δt = 2*calculate_correlation_time(;N = floor(Int, ρ * L^2), L = L)
        ρ = [0.2:0.2:4 ; [4.5, 5.0]]
        for ρ in ρ
            #The duration will be the number of data points times the time step plus 300 to eliminate the transient behaviour
            duration = number_of_data_points * Δt + 300
            N = floor(Int, ρ * L^2)
            config = @dict(L, N, duration, Δt, η)
            data = produce_or_load(execute_sim_dict, SimulationParameters(;config...), folder)[1]
            va = data["va"]
            @show mean(va[end-number_of_data_points:end])
            @show std(va[end-number_of_data_points:end])*2
            #We calculate the mean of the va for each noise value
            va = va[end-number_of_data_points:end]
            @show mean(va)
            @show std(va)*2 / (sqrt(length(va) - 1))
            @show std(va)*2 
            push!(va_mean, mean(va))
            push!(va_var, std(va, corrected = false))
            push!(va_error, std(va)*2 / sqrt(length(va) - 1)) 
        end
        #we save the data to csv
        safesave(datadir("va_vs_ρ", savename(@dict(L, η), "csv")), DataFrame(@strdict(ρ, va_mean, va_error, va_var)))
    end
end
