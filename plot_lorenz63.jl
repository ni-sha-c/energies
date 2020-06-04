include("test_lorenz63.jl")
using PyPlot
function plot_sensitivity()
	vsh, dJds = test_lss()
	fig, ax = subplots(1,1)
	vsh = vsh[:,:,3]'
	n = 0.005*(axes(vsh)[2] .- 1.)
	ax.plot(n, vsh[:,1], label=r"$v_{\rm sh,x}$")
	ax.plot(n, vsh[:,2], label=r"$v_{\rm sh,y}$")
	ax.plot(n, vsh[:,3], label=r"$v_{\rm sh,z}$")
	ax.set_xlabel("time")
	ax.set_ylabel("components of the shadowing direction")
	fig.legend()
end
