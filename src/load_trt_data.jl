struct TRTDataset
    full_data::DataFrame    # Full DataFrame with all data
    heating::DataFrame        # Heating phase segment
    cooling::DataFrame        # Cooling phase segment
end

"""
    load_trt_data(file_path; date_format=nothing, delim=',', header=1)

Imports experimental data from a thermal response test. The function assumes that the file is in a
file, with column as: [Time, Power [W], T_in [°C], T_out [°C]].
The `date_format` optional input allows to specify the format of the time column if it is in a
date-time format. If not provided, the function can auto-detect a supported date-time format
before falling back to numeric seconds.
# Arguments
    - `file_path`: Path to the data file
    - `date_format`: (Optional) Date format string for parsing time.
        - E.g., "yyyy-mm-dd HH:MM:SS"
    - `delim`: (Optional) Delimiter used in the CSV file (default is ',')
    - `header`: (Optional) Row number of the header (default is 1)
# Output
    - A DataFrame "trt" with columns: 
        - `:t`: Time array (either in seconds or DateTime)
        - `:power`: Heating power [W]
        - `:T_in`: Inlet temperature [°C]
        - `:T_out`: Outlet temperature [°C]
        - `:elapsed_time`: Elapsed time [s] (computed from DateTime when parsed, otherwise same as `:t`)
        - `:T_mean`: Mean fluid temperature[°C]
"""
function load_trt_data(file_path::String; date_format=nothing, delim=',', header::Int=1)
    # Read the CSV file into a DataFrame
    trt = CSV.read(file_path, DataFrame; delim=delim, header=header)
    
    # Rename for internal consistency
    rename!(trt, [1 => :t, 2 => :power, 3 => :T_in, 4 => :T_out])
    
    # Handle Date-Time parsing from explicit format or auto-detection.
    if !isnothing(date_format)
        trt.t = DateTime.(string.(trt.t), DateFormat(date_format))
        # Convert to elapsed seconds from start
        start_time = trt.t[1]
        trt.elapsed_time = map(t -> Dates.value(t - start_time) / 1000.0, trt.t)
    else
        trt.elapsed_time = Float64.(trt.t)  # Assume time is already in seconds
    end
    
    # Calculate Mean Fluid Temperature (Tf)
    trt.T_mean = mean_fluid_temperature(trt.T_in, trt.T_out, :arithmetic)
    return trt
end

"""
    decompose_trt(df, threshold=100.0)

Splits the TRT data into heating and cooling phases based on a power threshold.
# Arguments
    - `df`: DataFrame containing the TRT data with columns `:elapsed_time`, `:power`, and `:t_mean`
    - `threshold`: (Default 100.0) Threshold to distinguish between heating and cooling phases
"""
function decompose_trt(trt::DataFrame, threshold::Real=100.0)
    # Threshold-based heating identification.
    heating_mask = trt.power .> threshold
    heating_trt = trt[heating_mask, :]

    # Restitution/cooling is only considered after the heating phase.
    last_heating_idx = findlast(heating_mask)
    if isnothing(last_heating_idx)
        # If no heating phase is detected, return empty cooling dataset
        cooling_trt = trt[0:0, :]
    else
        # Cooling phase is after the last heating point with power below or equal to the threshold.
        after_heating = collect(1:nrow(trt)) .> last_heating_idx
        cooling_mask = (trt.power .<= threshold) .& after_heating
        cooling_trt = trt[cooling_mask, :]
    end
    
    return TRTDataset(trt, heating_trt, cooling_trt)
end