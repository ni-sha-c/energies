include("../examples/rijke.jl")
include("../src/clvs.jl")
using JLD
using OrdinaryDiffEq
function compute_les(expmt)
	d = N
	u = rand(d)
	nRunup = 1000000
	s = [7.0,0.2]
	u = Rijke_ODE(u, s, nRunup)
	nSteps = 10000
	dTu = zeros(d,d,nSteps)
	for n = 1:nSteps
		dTu[:,:,n] = dRijke(u, s, 1.e-8)
		u = Rijke_ODE(u, s, 1)
	end
	les, clv_trj = clvs(dTu, 20)
	save(string("../data/lyapunov/", 
		    "les_", string(expmt), 
		    ".jld"), "les", les)
end
function average_les()
	pmap(compute_les,897:1012)
end


