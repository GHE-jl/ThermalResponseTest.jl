using CairoMakie, DataFrames

function plot_trt_data(dataset::DataFrame)
    fig = Figure()
    
    # Plot Heating Phase
    ax1 = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Temperature (°C)")
    lines!(ax1, dataset.t, dataset.T_in, color=:red, linewidth=2, label="T_in")
    lines!(ax1, dataset.t, dataset.T_out, color=:blue, linewidth=2, label="T_out")
    lines!(ax1, dataset.t, dataset.T_mean, color=:green, linewidth=2, label="T_mean")
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Power (W)")
    lines!(ax2, dataset.t, dataset.power, color=:orange, linewidth=2, label="Power")
    axislegend(ax2, position=:rt)
    display(fig)
    
    return nothing
end