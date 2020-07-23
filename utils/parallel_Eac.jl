using Distributed
addprocs(16)
@everywhere include("plot_Eac.jl")
plot_acoustic_energy()
