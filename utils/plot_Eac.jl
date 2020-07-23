include("../examples/rijke.jl")
#using PyPlot
using OrdinaryDiffEq
using LinearAlgebra
using JLD
function plot_acoustic_energy()
n_samples = 500
Eac = zeros(n_samples)
beta = LinRange(5.0,9.0,n_samples)
nSteps = 20000
nRunup = 400000
first_LE = zeros(n_samples)
for i = 1:n_samples
	println("beta = ", beta[i])
	s = [beta[i], 0.2]
	u = Rijke_ODE(rand(N), s, nRunup)
	v = rand(N)
	v /= norm(v)
	for j = 1:nSteps	
		u = Rijke_ODE(u, s, 1)
		v = dRijke(u, s, 1.e-6)*v
		nv = norm(v)
		first_LE[i] += nv/nSteps
		v /= nv
		Eac[i] += norm(u[1:2*Ng])^2.0/4/nSteps
	end
	@show Eac[i]
end
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
save("../data/attractor/Eac.jld", "beta", beta, 
	 "Eac", Eac, "first_LE", first_LE)
return beta, Eac, first_LE

end
