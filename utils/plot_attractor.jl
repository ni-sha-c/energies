include("rijke.jl")
using PyPlot
n_runup = 40000
d = N
u0 = rand(d)
beta = 7.0
tau = 0.2
s = [beta, tau]
u = Rijke_ODE(u0, s, n_runup)[:,end]
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
ax.plot(qdot, uf, ".", ms=0.2)
ax.set_xlabel(L"\dot{q}",fontsize=18)
ax.set_ylabel(L"u_f",fontsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
ax.set_title(L"\beta = 7.0, \tau = 0.2", fontsize=18)
