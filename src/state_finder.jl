using Clustering

"""
    step_signal(x, steps)

Function that separated a vector signal "x" in n number of constant steps. This applies mainly
to help interprete noisy signal into constant values based on average abrupt changes.
Note: The k-means algorithm used to idenfity the changes can be unstable when too few steps are
used.
# Arguments
    - `x`: A vector
    - `steps`: The number of constant steps wanted in the output step-constant signal
# Output
    - A step-constant signal
"""
function step_signal(x::AbstractVector{<:Real}, steps::Integer)
    # Find group of same data using a k-mean cluster algorithm from Clustering.jl
    x_mat = reshape(x, 1, :)        # Convert to a 1 x n matrix
    res = kmeans(x_mat, steps)
    #result = kmedoids(x_mat, steps)

    # Compute mean for each cluster
    means = [sum(x[res.assignments .== k]) / length(x[res.assignments .== k]) for k in 1:steps]

    # Assign each value its cluster's mean
    return means[res.assignments]
end