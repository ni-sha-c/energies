using Distributed
addprocs(15)
@everywhere include("rijke_tangent_state_estimation.jl")
filenames = Array{String,1}(undef, nworkers())
for i = 1:nworkers()
	filenames[i] = string("../data/rijke_exp9_", string(i), ".jld") 
end
pmap(assimilate_parameter_and_trajectory, filenames)

