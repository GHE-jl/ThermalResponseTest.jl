"""
Script to test import options when loading TRT dataset
"""

includet("../src/load_trt_data.jl")
includet("../src/plot_trt_data.jl")

# Declare files to load (path from project root)
file1 = "trt_data/DataCL_TRT.csv"
file2 = "trt_data/DataCL.csv"
file3 = "trt_data/RAW_Mirabel1.txt"
file4 = "trt_data/TRT_IPL_1.txt"

trt1 = load_trt_data(file1; date_format=nothing, delim=',', header=1)
trt2 = load_trt_data(file2; date_format=nothing, delim=',', header=1)

# plot_trt_data(trt1)
# plot_trt_data(trt2)

Dataset1 = decompose_trt(trt1)
Dataset2 = decompose_trt(trt2)

plot_trt_data(Dataset1.heating)
plot_trt_data(Dataset2.heating)

plot_trt_data(Dataset1.cooling)
plot_trt_data(Dataset2.cooling)