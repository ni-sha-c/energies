using JLD
using PyPlot
function plot_sensitivity_cumsum()
X = load("../data/sensitivity_t2000_dt01_n500.jld")

dJds = X["dJdbeta"]
vsh = X["vsh"]
d, n, n_s = size(vsh)
t = 0.01*n
thres = 20
n_count = 0
mean_dJds = zeros(2)
var_dJds = zeros(2)
valid_samples = zeros(Int64, n_s)
for i = 1:n_s
	flag = 0
	for j = 1:n
		if norm(vsh[:,j,i]) > thres
			flag = 1
		end
	end
	if flag == 0
		valid_samples[i] = 1
		n_count += 1
		mean_dJds .+= dJds[:,i]
		var_dJds .+= dJds[:,i].*dJds[:,i]
	end
end

mean_dJds ./= n_count
var_dJds ./= n_count 
var_dJds .-= mean_dJds.*mean_dJds

fig, ax = subplots(1,2)
ax[1].xaxis.set_tick_params(labelsize=20)
ax[1].yaxis.set_tick_params(labelsize=20)
ax[2].xaxis.set_tick_params(labelsize=20)
ax[2].yaxis.set_tick_params(labelsize=20)

dJds = dJds[:,(!iszero).(valid_samples)]
n_s = n_count

ax[1].plot(t*(1:n_s), cumsum(dJds[1,:])./(1:n_s), "bo", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$")
ax[2].plot(t*(1:n_s), cumsum(dJds[2,:])./(1:n_s), "ro", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$")

ax[1].set_xlabel("Time",fontsize=20)
ax[1].set_ylabel("Shadowing sensitivities",fontsize=20)

ax[2].set_xlabel("Time",fontsize=20)
#ax[2].set_ylabel("Shadowing sensitivities",fontsize=20)
fig.legend(fontsize=20)
ax[1].grid(true)
ax[2].grid(true)


fig1, ax1 = subplots(1,1)
ax1.xaxis.set_tick_params(labelsize=20)
ax1.yaxis.set_tick_params(labelsize=20)

ax1.errorbar(1:n_count, dJds[1,:], sqrt(var_dJds[1]), marker="o", mfc="blue", mec="blue", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$")
ax1.set_xlabel("Sample #",fontsize=20)
ax1.set_ylabel("Shadowing sensitivities",fontsize=20)

fig1.legend(fontsize=20)
ax1.grid(true)

fig2, ax2 = subplots(1,1)
ax2.xaxis.set_tick_params(labelsize=20)
ax2.yaxis.set_tick_params(labelsize=20)

ax2.errorbar((1:n_count), dJds[2,:], sqrt(var_dJds[2]), marker="o", mfc="red", mec="red", color="red", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$")

ax2.set_xlabel("Sample #",fontsize=20)
ax2.set_ylabel("Shadowing sensitivities",fontsize=20)

fig2.legend(fontsize=20)
ax2.grid(true)
return fig1, ax1, fig2, ax2, dJds
end
