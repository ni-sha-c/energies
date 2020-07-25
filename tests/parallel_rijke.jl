using Distributed
addprocs(16)
@everywhere include("test_rijke.jl")
collect_sensitivities()
