include("../examples/rijke.jl")
using OrdinaryDiffEq
using PyPlot
using JLD
function plot_attractor(beta)
	n_runup = 400000
	d = N
	u0 = rand(d)
	tau = 0.2
	s = [beta, tau]
	u = Rijke_ODE(u0, s, n_runup)
	n_samples = 400000
	uf, qdot = zeros(n_samples), zeros(n_samples)
	u0 = copy(u)
	for k = 1:n_samples
		u .= Rijke_ODE(u0, s, 1)[:, end]
		uf[k] = dot(u[1:Ng], cjpixf)
		qdot[k] = beta*u[tNg]
		u0 .= u
	end
	save(string("../data/attractor/ufqdot_beta", beta, "_.jld"), 
		 "uf", uf, "qdot", qdot) 
	fig = figure(figsize=(10,6))
	ax = fig.add_subplot(111)
	ax.plot(qdot, uf, "b.", ms=0.2)
	ax.set_xlabel(L"\dot{q}",fontsize=46)
	ax.set_ylabel(L"u_f",fontsize=46)
	ax.xaxis.set_tick_params(labelsize=46)
	ax.yaxis.set_tick_params(labelsize=46)
	ax.set_title("\$\\beta = $(beta) \$", fontsize=46)
end
