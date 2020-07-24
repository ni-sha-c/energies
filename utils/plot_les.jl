include("../examples/rijke.jl")
include("../src/clvs.jl")
using JLD
using OrdinaryDiffEq
using PyPlot
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
function plot_les()
	n_files = 1012
	n_les = 20
	lyap_exps = zeros(n_les)
	for n = 1:n_files
		filename = string("../data/lyapunov/les_",
						  string(n), ".jld")
		X = load(filename)
		les = X["les"]
		lyap_exps += les/0.01/n_files
	end
	fig, ax = subplots(1,1)
	ax.plot(1:n_les,lyap_exps,"o")
	ax.set_xlabel("Index", fontsize=25)
	ax.set_ylabel("Lyapunov exponents", fontsize=25)
	ax.xaxis.set_tick_params(labelsize=25)
	ax.yaxis.set_tick_params(labelsize=25)
	ax.xaxis.set_ticks(1:n_les)
	grid(true)
	ax_inset = fig.add_axes([0.2,0.2,0.2,0.2])
	ax_inset.plot(1:3,lyap_exps[1:3],"o")
	ax_inset.xaxis.set_ticks(1:3)
	grid(true)
	ax_inset.xaxis.set_tick_params(labelsize=15)
	ax_inset.yaxis.set_tick_params(labelsize=15)
	return lyap_exps
end

