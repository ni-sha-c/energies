include("../examples/rijke.jl")
#include("../src/clvs.jl")
using PyPlot
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
function plot_rayleigh()
	filename = string("../data/attractor/",
					  "Jray.jld")
	X = load(filename)
	Jray = X["Jray"]
	beta = X["beta"]
	fig, ax = subplots(1,1)
	ax.plot(beta, Jray, "b.", ms=4.0)
	ax_in = fig.add_axes([0.6,0.6,0.2,0.2])
	ax_in.plot(beta[6.0 .<= beta .<= 7.2], 
			Jray[6.0 .<= beta .<= 7.2], "b.",ms=4.0)
	ax.xaxis.set_tick_params(labelsize=25)
	ax.yaxis.set_tick_params(labelsize=25)
	ax_in.xaxis.set_tick_params(labelsize=20)
	ax_in.yaxis.set_tick_params(labelsize=20)

	ax.set_xlabel(L"$\beta$", fontsize=25)
	ax.set_ylabel(L"$J_{\rm ray}$", fontsize=25)
	ax.grid(true)
	ax_in.grid(true)
end
