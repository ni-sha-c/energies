include("../examples/rijke.jl")
include("../src/clvs.jl")
using PyPlot
using OrdinaryDiffEq
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function compute_acoustic_energy()
    n_samples = 80
	d_u = 3
    Eac = SharedArray{Float64}(n_samples)
    les = SharedArray{Float64}(d_u, n_samples)
    les .= 0.
    Eac .= 0.
    beta = LinRange(5.0,9.0,n_samples)
    nSteps = 5000
    nRunup = 400000
    answer = @distributed for i = 1:n_samples
        println("beta = ", beta[i])
        s = [beta[i], 0.2]
        u = Rijke_ODE(rand(N), s, nRunup)
		du_trj = zeros(N, N, nSteps)
        for j = 1:nSteps	
        	u = Rijke_ODE(u, s, 1)
			du_trj[:,:,j] = dRijke(u, s, 1.e-6)
        	Eac[i] += norm(u[1:2*Ng])^2.0/4/nSteps
        end
		les[:,i], clv_trj = clvs(du_trj, d_u) 
        @show Eac[i]
    end
	wait(answer)
	save("../data/attractor/Eac.jld", "beta", beta, 
     "Eac", Eac, "les", les)
	return beta, Eac, les

end
function plot_acoustic_energy()
	X = load("../data/attractor/Eac.jld")
	Y = load("../data/attractor/Eac_and_beta.jld")
	Eac1 = X["Eac"]
	Eac2 = Y["Eac"]
	beta1 = X["beta"]
	beta2 = Y["beta"]
	first_LE = X["first_LE"]

	#=
	fig, ax = subplots(1,1)
    ax.plot(beta, Eac, "o")
    ax.xaxis.set_tick_params(labelsize=30)
    ax.yaxis.set_tick_params(labelsize=30)
    ax.set_xlabel(L"$\beta$",fontsize=30)
    ax.set_ylabel(L"$<E_{\rm ac}>$",fontsize=30)

    fig1, ax1 = subplots(1,1)
    ax1.plot(beta, first_LE, "o")
    ax1.xaxis.set_tick_params(labelsize=30)
    ax1.yaxis.set_tick_params(labelsize=30)
    ax1.set_xlabel(L"$\beta$",fontsize=30)
    ax1.set_ylabel(L"$\lambda^1$",fontsize=30)
	=#
end
