using Distributed
addprocs(16)
@everywhere include("plot_Eac.jl")
compute_acoustic_energy()
