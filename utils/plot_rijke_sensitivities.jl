using JLD
using PyPlot
function plot_tangent_sensitivities()
	n_files = 10000
	dJds_avg = zeros(2)
	dJds_var = zeros(2)
	dJds_arr = zeros(2,n_files)
	d = 28 
	nSteps = 2000
	vsh_arr = zeros(d, nSteps, n_files)
	for n = 1:480 
		filename = string("../data/",
						  "rijke_tangent_sensitivities/",
						  "vsh_and_dJds_",
						  string(n),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		#vsh = X["vsh"]
		dJds_arr[:,n] = dJds
		#vsh_arr[:,:,n] = vsh
	end
	for n = 481:1600
		filename = string("../data/",
						  "rijke_tangent_sensitivities/",
						  "vsh_and_dJds_",
						  string(n-480+1120),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		#vsh = X["vsh"]
		dJds_arr[:,n] = dJds
		#vsh_arr[:,:,n] = vsh
	end

	#indices = -10.0 .<= dJds_arr[1,:] .<= 10.0
	#n_files = sum(indices)
	#dJds_arr = dJds_arr[:,indices]
	X = load("../data/sensitivities/sensitivity_t2000_dt01_n500.jld")
	dJds = X["dJdbeta"]
	#indices = -10.0 .<= dJds[1,:] .<= 10.0
	#n_files += sum(indices)
	#dJds = dJds[:,indices]
	dJds_arr = [dJds_arr dJds]
	n_files = size(dJds_arr)[2]

	dJds_avg = sum(dJds_arr, dims=2)/n_files
	dJds_var = sum((dJds_arr .- dJds_avg).^2.0)/n_files 
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	
	ax.xaxis.set_tick_params(labelsize=28)
	ax1.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)
	
	ax.plot(1:n_files, dJds_arr[1,:], "b.", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$",ms=2.0)
	ax1.plot(1:n_files, dJds_arr[2,:], "r.", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$",ms=2.0)
	ax.plot(1:n_files, ones(n_files)*dJds_avg[1],"b",lw=2.0)
	ax1.plot(1:n_files, ones(n_files)*dJds_avg[2],"r",lw=2.0)
	ax.set_xlabel("Sample #",fontsize=28)
	ax.set_ylabel("Shadowing sensitivities",fontsize=28)
	ax.grid(true)
	ax1.set_xlabel("Sample #",fontsize=28)
	ax1.set_ylabel("Shadowing sensitivities",fontsize=28)
	ax1.grid(true)
	fig.legend(fontsize=28)
	return vsh_arr, dJds_arr, dJds_avg, dJds_var
end
function plot_tangent_sens_cum_avg()
	n_files = 9000
	dJds_avg = zeros(2)
	dJds_var = zeros(2)
	dJds_arr = zeros(2,n_files)
	nSteps = 2000
	#vsh_arr = zeros(d, nSteps, n_files)
	for n = 1:4000 
		filename = string("../data/",
						  "rijke_tangent_sensitivities/",
						  "withf/",
						  "vsh_and_dJds_",
						  string(n),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		#vsh = X["vsh"]
		dJds_arr[:,n] = dJds
		#vsh_arr[:,:,n] = vsh
	end
	for n = 4001:n_files 
		filename = string("../data/",
						  "rijke_tangent_sensitivities/",
						  "withf/",
						  "vsh_and_dJds_",
						  string(n + 1000),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		#vsh = X["vsh"]
		dJds_arr[:,n] = dJds
		#vsh_arr[:,:,n] = vsh
	end
	
	#indices = -10.0 .<= dJds_arr[1,:] .<= 10.0
	#n_files = sum(indices)
	#dJds_arr = dJds_arr[:,indices]
	X = load("../data/sensitivities/sensitivity_t2000_dt01_n500.jld")
	dJds = X["dJdbeta"]
	#indices = -10.0 .<= dJds[1,:] .<= 10.0
	#n_files += sum(indices)
	#dJds = dJds[:,indices]
	dJds_arr = [dJds_arr dJds]
	n_files = size(dJds_arr)[2]
	dJds_avg = sum(dJds_arr, dims=2)/n_files
	dJds_var = sum((dJds_arr .- dJds_avg).^2.0)/n_files 
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	
	ax.xaxis.set_tick_params(labelsize=28)
	ax1.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)
	T = 20.0 
	ax.plot(T*(1:n_files), 
			cumsum(dJds_arr[1,:])./(1:n_files),  
			"r.",  label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$",ms=4.0)
	#ax.set_ylim([-4.,-2.75])
	ax1.plot(T*(1:n_files), 
			 cumsum(dJds_arr[2,:])./(1:n_files), 
			 "b.", label=L"$\dfrac{d \langle J_{\rm ray}\rangle}{d \beta}$",ms=4.0)
	ax.plot(T*(1:n_files), ones(n_files)*dJds_avg[1],"r",lw=2.0)
	ax1.plot(T*(1:n_files), ones(n_files)*dJds_avg[2],"b",lw=2.0)
	ax.set_xlabel("Time",fontsize=28)
	ax.set_ylabel(L"$\dfrac{d\langle J_{\rm ac}\rangle}{d\beta}$",fontsize=28)
	ax.grid(true)
	ax1.set_xlabel("Time",fontsize=28)
	ax1.set_ylabel(L"$\dfrac{d\langle J_{\rm ray}\rangle}{d\beta}$",fontsize=28)
	#ax1.set_ylim([0,0.05])
	ax1.grid(true)
	lgnd = fig.legend(loc=(0.9,0.7),fontsize=28)
	lgnd.legendHandles[1]._legmarker.set_markersize(10)
	lgnd = fig1.legend(loc=(0.9,0.2),fontsize=28)
	lgnd.legendHandles[1]._legmarker.set_markersize(10)
	
	return dJds_arr
end

function plot_adjoint_sensitivities()
	n_files = 10000
	dJds_avg = zeros(2)
	dJds_var = zeros(2)
	dJds_arr = zeros(2,n_files)
	d = 28 
	nSteps = 2800
	vsh_arr = zeros(d, nSteps, n_files)
	for n = 1:n_files 
		filename = string("../data/",
						  "rijke_adjoint_sensitivities/",
						  "vsh_and_dJds_",
						  string(n),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		vsh = X["vsh"]
		dJds *= nSteps
		dJds_arr[:,n] = dJds
		vsh_arr[:,:,n] = vsh
	end
	indices = -10.0 .<= dJds_arr[1,:] .<= 10.0
	n_files = sum(indices)
	dJds_arr = dJds_arr[:, indices]
	dJds_avg = sum(dJds_arr, dims=2)/n_files
	dJds_var = sum((dJds_arr .- dJds_avg).^2.0, dims=2)/n_files 
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax1.xaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)
	
	ax.plot(1:n_files, dJds_arr[1,:], "b.", label=L"$\dfrac{d \langle E_{\rm ac}\rangle}{d \beta}$",ms=2.0)
	ax1.plot(1:n_files, dJds_arr[2,:], "r.", label=L"$\dfrac{d \langle E_{\rm ac}\rangle}{d \tau}$",ms=2.0)
	ax.plot(1:n_files, ones(n_files)*dJds_avg[1],"b",lw=2.0)
	ax1.plot(1:n_files, ones(n_files)*dJds_avg[2],"r",lw=2.0)
	ax.set_xlabel("Sample #",fontsize=28)
	ax.set_ylabel("Shadowing sensitivities",fontsize=28)
	ax.grid(true)
	ax1.set_xlabel("Sample #",fontsize=28)
	ax1.set_ylabel("Shadowing sensitivities",fontsize=28)
	ax1.grid(true)
	fig.legend(fontsize=28)
	return dJds_arr
end
function plot_adjoint_sens_cum_avg()
	n_files = 10000
	dJds_avg = zeros(2)
	dJds_var = zeros(2)
	dJds_arr = zeros(2,n_files)
	nSteps = 2000
	#vsh_arr = zeros(d, nSteps, n_files)
	for n = 1:n_files 
		filename = string("../data/",
						  "rijke_adjoint_sensitivities/",
						  "6.8new/",
						  "vsh_and_dJds_",
						  string(n),
						  ".jld")

		X = load(filename)
		dJds = X["dJds"]
		#vsh = X["vsh"]
		dJds *= nSteps
		dJds_arr[:,n] = dJds
		#vsh_arr[:,:,n] = vsh
	end
	#indices = -20.0 .<= dJds_arr[1,:] .<= 20.0
	#n_files = sum(indices)
	#dJds_arr = dJds_arr[:, indices]
	dJds_avg = sum(dJds_arr, dims=2)/n_files
	dJds_var = sum((dJds_arr .- dJds_avg).^2.0, dims=2)/n_files 
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax1.xaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)
	T = 20.
	ax.plot(T*(1:n_files), 
			cumsum(dJds_arr[1,:])./(1:n_files), "r.", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \beta}$",ms=4.0)
	ax1.plot(T*(1:n_files), 
			 cumsum(dJds_arr[2,:])./(1:n_files), "g.", label=L"$\dfrac{d \langle J_{\rm ac}\rangle}{d \tau}$",ms=4.0)
	ax.plot(T*(1:n_files), ones(n_files)*dJds_avg[1],"r",lw=2.0)
	ax1.plot(T*(1:n_files), ones(n_files)*dJds_avg[2],"g",lw=2.0)

	ax.set_ylim([-4.9,-3.1])
	ax1.set_ylim([100,200])
	ax.set_xlabel("Time",fontsize=28)
	ax.set_ylabel(L"$\dfrac{d\langle J_{\rm ac}\rangle}{d\beta}$",fontsize=28)
	ax.grid(true)
	ax1.set_xlabel("Time",fontsize=28)
	ax1.set_ylabel(L"$\dfrac{d\langle J_{\rm ac}\rangle}{d\tau}$",fontsize=28)
	#ax1.set_ylabel("Shadowing sensitivities",fontsize=28)
	ax1.grid(true)
	#lgnd = fig.legend(loc=(0.5,0.6),fontsize=28)
	#lgnd.legendHandles[1]._legmarker.set_markersize(10)
	#lgnd = fig1.legend(loc=(0.5,0.6),fontsize=28)
	#lgnd.legendHandles[1]._legmarker.set_markersize(10)
	plt.xticks(LinRange(0,200000,5))
	#plt.locator_params(axis='x', nbins=6)

	return dJds_arr
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
