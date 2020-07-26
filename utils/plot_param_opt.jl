using JLD
using PyPlot
using Distributed
using SharedArrays
function plot_optimization_path()
		filename = string("../data/",
						  "optimization/",
						  "beta_dJds",
						  ".jld")

		filename1 = string("../data/",
						  "optimization/",
						  "beta_dJds_1",
						  ".jld")



		X = load(filename)
		dJds = X["dJds_path"]
		beta_path = X["beta_path"]
		
		Y = load(filename1)
		dJds1 = Y["dJds_path"]
		beta1_path = Y["beta_path"]

		X = load("../data/optimization/Eac_along_path.jld")
		beta_path_r = X["beta"]
		beta1_path_r = X["beta1"]
		Eac_path = X["Eac"]
		Eac1_path = X["Eac1"]

		return beta_path_r, beta1_path_r, Eac_path, Eac1_path
		#=
		X = load("../data/attractor/more_les/Eac.jld")
		Y = load("../data/attractor/more_les/Eac1.jld")
		Z = load("../data/attractor/Eac.jld")
		W = load("../data/attractor/Eac_and_beta.jld")
		Eac1 = X["Eac"]
		Eac2 = Y["Eac"]
		beta1 = X["beta"]
		beta2 = Y["beta"]
		Eac3 = Z["Eac"]
		Eac4 = W["Eac"]
		beta3 = Z["beta"]
		beta4 = W["beta"]
		
		ll = 5.5
		ul = 7.4

		Eac1 = Eac1[ll .<= beta1 .<= ul]
		Eac2 = Eac2[ll .<= beta2 .<= ul]
		Eac3 = Eac3[ll .<= beta3 .<= ul]
		Eac4 = Eac4[ll .<= beta4 .<= ul]
	
		beta1 = beta1[ll .<= beta1 .<= ul]
		beta2 = beta2[ll .<= beta2 .<= ul]
		beta3 = beta3[ll .<= beta3 .<= ul]
		beta4 = beta4[ll .<= beta4 .<= ul]



		fig, ax = subplots(1,1)
		ax.xaxis.set_tick_params(labelsize=25)
		ax.yaxis.set_tick_params(labelsize=25)
	
		ax.plot(beta1, Eac1, "r.",ms=2.0)
		ax.plot(beta2, Eac2, "r.",ms=2.0)
		ax.plot(beta3, Eac3, "r.",ms=2.0)
		ax.plot(beta4, Eac4, "r.",ms=2.0)
		#=
		ax.plot(1:n_files, ones(n_files)*dJds_avg[1],"b",lw=2.0)
	ax.set_xlabel("Sample #",fontsize=25)
	ax.set_ylabel("Shadowing sensitivities",fontsize=25)
	=#
	ax.grid(true)
	
	return beta, dJds
	=#
end
