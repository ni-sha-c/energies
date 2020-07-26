using Distributed
addprocs(16)
@everywhere include("optimize.jl")
#compute_acoustic_energy()
compute_Eac_along_path()
