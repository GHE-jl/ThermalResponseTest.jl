"""
    TRTDataset(data, heating, recovery)

Data segments of a **conventional** thermal response test: full dataset (`data`), heating phase
(`heating`), and recovery phase (`recovery`), as produced by [`decompose_trt`](@ref).

An unconventional TRT (e.g. several on/off cycles) does not fit this single heating/recovery shape;
interpret it with the vector-based (`t`, `T`, `q`, ...) methods directly on `load_trt_data`'s
DataFrame instead.
"""
struct TRTDataset
    data::DataFrame
    heating::DataFrame
    recovery::DataFrame
end

"""
    load_trt_data(file_path; date_format=nothing, delim=',', header=1, dt_tol=0.1,
        trim_recirculation=true, threshold=100.0)

Imports experimental data from a thermal response test: `[Time, Power [W], T_in [°C], T_out [°C]]`.
`date_format` parses a date-time column (e.g. `"yyyy-mm-dd HH:MM:SS"`); otherwise time is assumed
already numeric (seconds).

The output starts at `t = 0` with `power = 0`, giving the undisturbed ground temperature `T0` before
heating is turned on at the next sample. Non-uniform time steps are PCHIP-interpolated onto a 
uniform grid if they drift by more than `dt_tol` from the median. Any recirculation phase before
heating is dropped, keeping only the sample right before heating turns on, with `:elapsed_time` 
rebased so that sample is `t = 0`.
# Arguments
    - `file_path`: Path to the data file
    - `date_format`, `delim`, `header`: CSV parsing options
    - `dt_tol`: Max deviation (s) from the median time step before interpolating
    - `trim_recirculation`: Drop any pre-heating recirculation phase (default `true`)
    - `threshold`: Heating power threshold used to find the start of heating (W)
# Output
    - DataFrame with `:elapsed_time`, `:power`, `:T_in`, `:T_out`, `:T_mean`
"""
function load_trt_data(file_path::String; date_format=nothing, delim=',', header::Int=1,
    dt_tol::Real=0.1, trim_recirculation::Bool=true, threshold::Real=100.0)

    # Load the data
    trt = CSV.read(file_path, DataFrame; delim=delim, header=header)
    rename!(trt, [1 => :t, 2 => :power, 3 => :T_in, 4 => :T_out])

    # Check for DateTime format
    if !isnothing(date_format)
        trt.t = DateTime.(string.(trt.t), DateFormat(date_format))
        start_time = trt.t[1]
        trt.elapsed_time = map(t -> Dates.value(t - start_time) / 1000.0, trt.t)
    else
        trt.elapsed_time = Float64.(trt.t)
    end

    # Uniform time steps.
    steps = diff(trt.elapsed_time)
    dt = sort(steps)[cld(length(steps), 2)]
    if maximum(abs.(steps .- dt)) > dt_tol
        @warn "Non-uniform time steps, interpolating onto a uniform grid (dt = $dt s)."
        t_uniform = collect(trt.elapsed_time[1]:dt:trt.elapsed_time[end])
        trt = DataFrame(
            elapsed_time = t_uniform,
            power = pchip_interpolation(trt.elapsed_time, trt.power, t_uniform),
            T_in  = pchip_interpolation(trt.elapsed_time, trt.T_in, t_uniform),
            T_out = pchip_interpolation(trt.elapsed_time, trt.T_out, t_uniform),
        )
    end

    # Recirculation trim: keep only the sample right before heating turns on, and rebase to t = 0.
    if trim_recirculation
        first_heating_idx = findfirst(>(threshold), trt.power)
        if first_heating_idx == 1
            @warn "No pre-heating sample found; T0 cannot be determined from this file."
        elseif !isnothing(first_heating_idx)
            trt = trt[(first_heating_idx - 1):end, :]
        end
        trt.elapsed_time = trt.elapsed_time .- trt.elapsed_time[1]
        trt.elapsed_time[1] == 0 ? trt.elapsed_time[1] = 1.0 : nothing
    end

    trt.T_mean = mean_fluid_temperature(trt.T_in, trt.T_out, :arithmetic)
    return trt
end

"""
    decompose_trt(trt, threshold=100.0)

Splits a **conventional** TRT (one heating phase, then one recovery phase) into a
[`TRTDataset`](@ref) based on a power threshold.
# Arguments
    - `trt`: DataFrame from `load_trt_data`
    - `threshold`: Heating power threshold separating heating from recovery
# Output
    - `TRTDataset` (`heating`/`recovery` also carry `:t_rel`, time since the start of heating)
"""
function decompose_trt(trt::DataFrame, threshold::Real=100.0)
    heating_mask = trt.power .> threshold
    heating_trt = copy(trt[heating_mask, :])

    last_heating_idx = findlast(heating_mask)
    if isnothing(last_heating_idx)
        recovery_trt = copy(trt[0:0, :])
    else
        after_heating = collect(1:nrow(trt)) .> last_heating_idx
        recovery_mask = (trt.power .<= threshold) .& after_heating
        recovery_trt = copy(trt[recovery_mask, :])
    end

    # Rebase both phases to t = 0 at the moment heating starts.
    if nrow(heating_trt) > 0 && nrow(trt) > 1
        Δt = trt.elapsed_time[2] - trt.elapsed_time[1]
        t0 = heating_trt.elapsed_time[1] - Δt
    else
        t0 = 0.0
    end
    heating_trt.t_rel = heating_trt.elapsed_time .- t0
    recovery_trt.t_rel = recovery_trt.elapsed_time .- t0

    data = copy(trt)
    nrow(data) > 0 && data.elapsed_time[1] == 0 && (data.elapsed_time[1] = 1.0)

    return TRTDataset(data, heating_trt, recovery_trt)
end
