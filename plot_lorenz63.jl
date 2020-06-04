include("test_lorenz63.jl")
using PyPlot
function plot_sensitivity()
	vsh, dJds = test_lss()
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	vsh = vsh[:,100:end,10]'
	n = 0.005*(axes(vsh)[1] .- 1.)
	n_samples = size(dJds)[1]
	dJds_avg = sum(dJds)/n_samples
	var_dJds = sum(x -> x^2, dJds .- dJds_avg)/n_samples
	ax.plot(n, vsh[:,1], label=L"v_{\rm sh,x}")
	ax.plot(n, vsh[:,2], label=L"v_{\rm sh,y}")
	ax.plot(n, vsh[:,3], label=L"v_{\rm sh,z}")
	ax.set_xlabel("time",fontsize=30)
	ax.set_ylabel("components of the shadowing direction",
					  fontsize=18)
	ax.xaxis.set_tick_params(labelsize=18)
	ax.yaxis.set_tick_params(labelsize=18)
	fig.legend()
	ax1.set_ylabel("sensitivities from NILSS",fontsize=18)
	ax1.set_xlabel("sample number",fontsize=18)
	ax1.set_title("trajectory length = 5000", fontsize=18)
	ax1.xaxis.set_tick_params(labelsize=18)
	ax1.yaxis.set_tick_params(labelsize=18)
	ax1.grid(true)
	ax1.errorbar(x=1:n_samples,y=dJds,yerr=sqrt(var_dJds),
			 linestyle="none",ms=4)
	ax1.plot(sum(dJds)/n_samples*ones(n_samples),"--")
	#mean value: 0.905328466, variance = 0.042
end
function plot_condition_number()
	dJds, condnum = test_condition_number()
	fig, ax = subplots(1,1)
	n_samples = size(dJds)[1]
	ax.plot(1:n_samples, condnum, ".", ms=4.0)
	ax.set_xlabel("trajectory length",fontsize=18)
	ax.set_ylabel("condition number",
					  fontsize=18)
	ax.xaxis.set_tick_params(labelsize=18)
	ax.yaxis.set_tick_params(labelsize=18)
end
