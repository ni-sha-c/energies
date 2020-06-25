include("../examples/rijke.jl")
using PyPlot
function plot_acoustic_energy()
N = 30
n_samples = 150
Eac = zeros(n_samples)
beta = LinRange(5.6,7.6,n_samples)
nSteps = 20000
nRunup = 60000
for i = 1:n_samples
	println("beta = ", beta[i])
	s = [beta[i], 0.2]
	u = Rijke_ODE(rand(N), s, nRunup)
	for j = 1:nSteps	
		u = Rijke_ODE(u, s, 1)
		Eac[i] += norm(u[1:2*Ng])^2.0/4/nSteps
	end
	@show Eac[i]
end
fig, ax = subplots(1,1)
ax.plot(beta, Eac, "o")
ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.set_xlabel(L"$\beta$",fontsize=20)
ax.set_ylabel(L"$<E_{\rm ac}>$",fontsize=20)
return fig, ax, beta, Eac
end
