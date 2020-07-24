using JLD
using PyPlot
function plot_adjoint_sensitivities()
	n_files = 1120
	dJds_avg = zeros(2)
	dJds_var = zeros(2)
	dJds_arr = zeros(2,n_files)
	d = 30 
	n = 2000
	vsh_arr = zeros(d, n, n_files)
	for n = 1:n_files 
		filename = string("../data/",
						  "rijke_adjoint_sensitivities/",
						  "vsh_and_dJds_",
						  string(n),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		vsh = X["vsh"]
		dJds_arr[:,n] = dJds
		dJds_avg += dJds/n_files
		dJds_var += dJds.*dJds/n_files
		vsh_arr[:,:,n] = vsh
	end
	dJds_var .-= dJds_avg.*dJds_avg

	fig, ax = subplots(1,2)
	ax[1].xaxis.set_tick_params(labelsize=25)
	ax[1].yaxis.set_tick_params(labelsize=25)
	ax[2].xaxis.set_tick_params(labelsize=25)
	ax[2].yaxis.set_tick_params(labelsize=25)
	
	ax[1].plot(1:n_files, dJds_arr[1,:], "b.", label=L"$\dfrac{d \langle E_{\rm ac}\rangle}{d \beta}$")
	ax[2].plot(1:n_files, dJds_arr[2,:], "r.", label=L"$\dfrac{d \langle E_{\rm ac}\rangle}{d \tau}$")
	ax[1].plot(1:n_files, ones(n_files)*dJds_avg[1],"b--")
	ax[2].plot(1:n_files, ones(n_files)*dJds_avg[2],"r--")
	ax[1].set_xlabel("Sample #",fontsize=25)
	ax[1].set_ylabel("Shadowing sensitivities",fontsize=25)
	ax[1].grid(true)
	ax[2].set_xlabel("Sample #",fontsize=25)
	ax[2].set_ylabel("Shadowing sensitivities",fontsize=25)
	ax[2].grid(true)
	fig.legend(fontsize=25)
	return vsh_arr, dJds_arr, dJds_avg, dJds_var
end

#=
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
=#
