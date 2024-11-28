using NearestNeighbors
using BenchmarkTools

data = rand(3, 10^4)
k = 3
point = @SVector rand(3)

kdtree = KDTree(data)
idxs, dists = knn(kdtree, point, k, true)

@benchmark nn($kdtree, $point)

# idxs
# # 3-element Array{Int64,1}:
# #  4683
# #  6119
# #  3278

# dists
# # 3-element Array{Float64,1}:
# #  0.039032201026256215
# #  0.04134193711411951
# #  0.042974090446474184

# # Multiple points
# points = rand(3, 4)
# idxs, dists = knn(kdtree, points, k, true)

# idxs
# # 4-element Array{Array{Int64,1},1}:
# #  [3330, 4072, 2696]
# #  [1825, 7799, 8358]
# #  [3497, 2169, 3737]
# #  [1845, 9796, 2908]

# # dists
# # 4-element Array{Array{Float64,1},1}:
# #  [0.0298932, 0.0327349, 0.0365979]
# #  [0.0348751, 0.0498355, 0.0506802]
# #  [0.0318547, 0.037291, 0.0421208]
# #  [0.03321, 0.0360935, 0.0411951]

# # Static vectors
# using StaticArrays
# v = @SVector[0.5, 0.3, 0.2];

# idxs, dists = knn(kdtree, v, k, true)

# idxs
# # 3-element Array{Int64,1}:
# #   842
# #  3075
# #  3046

# dists
# # 3-element Array{Float64,1}:
# #  0.04178677766255837
# #  0.04556078331418939
# #  0.049967238112417205

# # Preallocating input results
# idxs, dists = zeros(Int32, k), zeros(Float32, k)
# knn!(idxs, dists, kdtree, v, k)