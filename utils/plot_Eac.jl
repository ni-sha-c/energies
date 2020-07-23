include("../examples/rijke.jl")
#using PyPlot
using OrdinaryDiffEq
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function plot_acoustic_energy()
n_samples = 200
Eac = SharedArray{Float64}(n_samples)
first_LE = SharedArray{Float64}(n_samples)
first_LE .= 0.
Eac .= 0.
beta = LinRange(2.0,5.0,n_samples)
nSteps = 10000
nRunup = 1000000
answer = @distributed for i = 1:n_samples
	println("beta = ", beta[i])
	s = [beta[i], 0.2]
	u = Rijke_ODE(rand(N), s, nRunup)
	v = rand(N)
	v /= norm(v)
	for j = 1:nSteps	
		u = Rijke_ODE(u, s, 1)
		v = dRijke(u, s, 1.e-6)*v
		nv = norm(v)
		first_LE[i] += log(nv)/nSteps
		v /= nv
		Eac[i] += norm(u[1:2*Ng])^2.0/4/nSteps
	end
	@show Eac[i]
end
wait(answer)
save("../data/attractor/Eac_2_to_5.jld", "beta", beta, 
	 "Eac", Eac, "first_LE", first_LE)
return beta, Eac, first_LE

end
