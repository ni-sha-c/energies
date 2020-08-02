include("../examples/rijke.jl")
using OrdinaryDiffEq
using PyPlot
n_runup = 100000
d = N
u0 = rand(d)
beta = 0.7
tau = 0.2
s = [beta, tau]
u = Rijke_ODE(u0, s, n_runup)
n_samples = 1000000
uf, qdot = zeros(n_samples), zeros(n_samples)
u0 = copy(u)
for k = 1:n_samples
	u .= Rijke_ODE(u0, s, 1)[:, end]
	uf[k] = dot(u[1:Ng], cjpixf)
	qdot[k] = beta*u[tNg]
	u0 .= u
end
fig, ax = subplots(1,1)
ax.plot(qdot, uf, "g.", ms=0.2)
ax.set_xlabel(L"\dot{q}",fontsize=28)
ax.set_ylabel(L"u_f",fontsize=28)
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.set_title(L"\beta = 8.0, \tau = 0.2", fontsize=28)
