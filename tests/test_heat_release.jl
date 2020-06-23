include("rijke.jl")
using PyPlot
N = 18
u0Init = zeros(N)
u0Init[1] = 1.
s = [7.0, 0.2]
nRunUp = 40000
#u0 = Rijke(u0Init, s, nRunUp)
t = LinRange(-1.01,-0.99,500)
heat_sqrt = qfun_sqrt.(t)
heat_poly = qfun_poly.(t)
fig, ax = subplots(1,1)
ax.plot(t, heat_sqrt, label=L"$\dot{q}$, using sqrt", ".")
ax.plot(t, heat_poly, label=L"$\dot{q}$, using poly", ".")
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
ax.set_xlabel("delayed velocity",fontsize=18)
fig.legend(fontsize=18)


