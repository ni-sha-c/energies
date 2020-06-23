include("rijke.jl")
N = 30
n_samples = 50
Eac = zeros(n_samples)
beta = LinRange(5.0,9.0,n_samples)
u = zeros(N)
u[1] = 1.
nSteps = 500000
for i = 1:n_samples
		Eac[i] = rijke(u, beta[i], nSteps)
end
fig, ax = subplots(1,1)
ax.plot(beta, Eac, 'o')
ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.set_xlabel(L"\beta",fontsize=20)
ax.set_ylabel(L"<E_{\rm ac}>",fontsize=20)

