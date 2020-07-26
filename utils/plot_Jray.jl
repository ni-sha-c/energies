include("../examples/rijke.jl")
#include("../src/clvs.jl")
#using PyPlot
using OrdinaryDiffEq
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function compute_rayleigh_criterion(n_samples)
    d_u = 3
    Jray = SharedArray{Float64}(n_samples)
    Jray .= 0.
    beta = LinRange(6.0,9.0,n_samples)
    nSteps = 20000
    nRunup = 400000
    answer = @distributed for i = 1:n_samples
    	println("beta = ", beta[i])
    	s = [beta[i], 0.2]
    	u = Rijke_ODE(rand(N), s, nRunup)
    	for j = 1:nSteps    
        	u = Rijke_ODE(u, s, 1)
        	Jray[i] += sum(u[1:Ng].*cjpixf)/nSteps
    	end
    	@show Jray[i]
    end
    wait(answer)
    save("../data/attractor/Jray.jld", "beta", beta, 
     "Jray", Jray)
end

