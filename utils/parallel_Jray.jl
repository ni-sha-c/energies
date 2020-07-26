using Distributed
addprocs(16)
@everywhere include("plot_Jray.jl")
compute_rayleigh_criterion(1120)
