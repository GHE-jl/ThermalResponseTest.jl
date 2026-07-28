module ThermalResponseTest

# GroundHeatExchanger re-exports the GroundResponse.jl ground models (ils, ics, fls, mils, mfls)
# and the temporal-superposition `convolution`, which the model inversions are built on.
using GroundHeatExchanger

using CSV, DataFrames, Dates
using LinearAlgebra
using Clustering
using Optimization, OptimizationOptimJL
using FiniteDiff   # loads the AutoFiniteDiff backend used by the inversions

# Source files (included once, in dependency order)
include("utils.jl")                     # mean fluid temperature, derivative, time helpers
include("trt_data.jl")                  # TRTDataset, load_trt_data, decompose_trt
include("first_order_approximation.jl") # ILS first-order-approximation interpretation methods
include("trt_inversion.jl")             # Optimization.jl model inversions
include("state_finder.jl")              # step_signal (variable-power segmentation)

# Data loading and pre-processing
export TRTDataset, load_trt_data, decompose_trt

# First-order approximation (ILS) interpretation — Pasquier (2018)
export fit_ils_foa_T            # UFOA-T-H : temperature, heating phase (Eq. 3)
export fit_ils_foa_T_recovery  # UFOA-T-R : temperature, recovery phase (Eq. 13)
export fit_ils_foa_dT          # CFOA-Ṫ-H : temperature derivative, heating phase (Eqs. 8–10)
export fit_ils_foa_dT_recovery # CFOA-Ṫ-R : temperature derivative, recovery phase (Eqs. 14–18)

# Model inversion (Optimization.jl + Optim.jl), ground conductivity as main unknown
export fit_ground_response                       # generic core (any AbstractGroundModel)
export fit_ils, fit_ics, fit_fls, fit_mils, fit_mfls

# Variable-power segmentation
export step_signal

# Utilities
export mean_fluid_temperature, centered_finite_difference, bourdet_derivative, critical_time,
    residence_time

end
