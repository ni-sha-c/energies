using Distributed
addprocs(16)
@everywhere include("plot_les.jl")
average_les()
