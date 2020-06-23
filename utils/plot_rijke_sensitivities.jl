using JLD
using PyPlot
X = load("rijke_shadowing_sens.jld")

n_arr = X["n"]
n = length(n_arr)

dJds = X["dJds"][:,1:n]
var_dJds = X["var_dJds"][:,1:n]

fig, ax = subplots(1,1)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)

ax.plot(n_arr*0.01, dJds[1,:], "bo", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$")
ax.plot(n_arr*0.01, dJds[2,:], "ro", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$")

ax.set_xlabel("time segment length",fontsize=18)
ax.set_ylabel("Shadowing sensitivities",fontsize=18)

fig.legend(fontsize=18)
ax.grid(true)

fig1, ax1 = subplots(1,1)
ax1.xaxis.set_tick_params(labelsize=18)
ax1.yaxis.set_tick_params(labelsize=18)

ax1.errorbar(n_arr*0.01, dJds[1,:], sqrt.(var_dJds[1,:]), marker="o", mfc="blue", mec="blue", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$")
ax1.set_xlabel("time segment length",fontsize=18)
ax1.set_ylabel("Shadowing sensitivities",fontsize=18)

fig1.legend(fontsize=18)
ax1.grid(true)

fig2, ax2 = subplots(1,1)
ax2.xaxis.set_tick_params(labelsize=18)
ax2.yaxis.set_tick_params(labelsize=18)

ax2.errorbar(n_arr*0.01, dJds[2,:], sqrt.(var_dJds[2,:]), marker="o", mfc="red", mec="red", color="red", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$")

ax2.set_xlabel("time segment length",fontsize=18)
ax2.set_ylabel("Shadowing sensitivities",fontsize=18)

fig2.legend(fontsize=18)
ax2.grid(true)

