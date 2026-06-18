"""
Script to demonstrate the application of 
"""

using DelimitedFiles
using Clustering
using CairoMakie

includet("../src/Convolutions.jl")

# Load the heating load in the reference file
filename = "tests/TRT_StepTest1_QVar.csv"
DATA = readdlm(filename, ',', skipstart=1)
Q = vec(DATA[:, 2]) # Heating load in W

# Test clustering
function test_clustering(Q::Vector{Float64}, n_clusters::Int)
    # Apply k-means clustering to the heating load data
    Q = reshape(Q, 1, :)  # Reshape to a 1 x n matrix
    kmeans_result = kmeans(Q, n_clusters; maxiter=100, display=:none)
    
    # Extract cluster centers and assignments
    cluster_centers = kmeans_result.centers
    assignments = kmeans_result.assignments
    
    return cluster_centers, assignments
end

# Test step_signal function
@time x = step_signal(Q, 100)

f = Figure(; size = (16 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1]; xlabel = L"$t$ (s)", ylabel = L"$Q$ (W)")
lines!(ax, 1:length(Q), Q; color = :blue, label = "Original signal")
lines!(ax, 1:length(x), x; color = :red, label = "Step signal (100 steps)")
axislegend(ax; position = :rt)
display(f)

# Test the clustering function
n_clusters = 100
centers, assignments = test_clustering(Q, n_clusters)