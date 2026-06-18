module ThermalResponseTest

using Revise
include("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")    # Load the GroundHeatExchanger
using .GroundHeatExchanger

# Load files
includet("load_trt_data.jl")
includet("first_order_approximation.jl")
includet("trt_inversion.jl")
includet("plot_trt_data.jl")
includet("state_finder.jl")
includet("utils.jl")

export TRTData, TRTDataset, load_trt_data, decompose_trt

export fit_ils_foa_T, fit_ils_foa_dT

export fit_trt_parameters

export plot_trt_data

export step_signal

export mean_fluid_temperature, centered_finite_difference, critical_time,
    residence_time, residence_time_indice

end