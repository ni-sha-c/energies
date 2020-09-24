using JLD
using PyPlot
using Distributed
using SharedArrays
function plot_optimization_path()
		#filename = string("../data/",
		#				  "optimization/",
		#				  "beta_dJds",
		#				  ".jld")

		#filename1 = string("../data/",
		#				  "optimization/",
		#				  "beta_dJds_1",
		#				  ".jld")



		#X = load(filename)
		#dJds = X["dJds_path"]
		#beta_path = X["beta_path"]
		
		#Y = load(filename1)
		#dJds1 = Y["dJds_path"]
		#beta1_path = Y["beta_path"]

		X = load("../data/optimization/Eac_along_path_05.jld")
		beta = X["beta"]
		#beta1 = X["beta1"]
		Eac = X["Eac"]
		#Eac1 = X["Eac1"]

		X = load("../data/attractor/more_les/Eac.jld")
		Y = load("../data/attractor/more_les/Eac1.jld")
		Z = load("../data/attractor/Eac.jld")
		W = load("../data/attractor/Eac_and_beta.jld")
		J1 = X["Eac"]
		J2 = Y["Eac"]
		s1 = X["beta"]
		s2 = Y["beta"]
		J3 = Z["Eac"]
		J4 = W["Eac"]
		s3 = Z["beta"]
		s4 = W["beta"]
		
		ll = 4.5
		ul = 7.5

		J1 = J1[ll .<= s1 .<= ul]
		J2 = J2[ll .<= s2 .<= ul]
		J3 = J3[ll .<= s3 .<= ul]
		J4 = J4[ll .<= s4 .<= ul]
	
		s1 = s1[ll .<= s1 .<= ul]
		s2 = s2[ll .<= s2 .<= ul]
		s3 = s3[ll .<= s3 .<= ul]
		s4 = s4[ll .<= s4 .<= ul]



		fig, ax = subplots(1,1)
		ax.xaxis.set_tick_params(labelsize=36)
		ax.yaxis.set_tick_params(labelsize=36)
	
		ax.plot(s1, J1, "r.",ms=6.0)
		ax.plot(s2, J2, "r.",ms=6.0)
		ax.plot(s3, J3, "r.",ms=6.0)
		ax.plot(s4, J4, "r.",ms=6.0)
		
		ax.plot(beta[1:10], Eac[1:10], "bP", ms=8.0)
		#ax.plot(beta1[1:6], Eac1[1:6], "kP", ms=8.0)
		#ax.plot(6.5, 23.1, "bP", ms=8.0)
		
		ax.grid(true)
		ax.set_xlabel(L"$\beta$", fontsize=36)
		ax.set_ylabel(L"$\langle J_{\rm ac}\rangle$", 
					  rotation=0,fontsize=36)


		ll = 6.0
		ul = 7.3

		J1 = J1[ll .<= s1 .<= ul]
		J2 = J2[ll .<= s2 .<= ul]
		J3 = J3[ll .<= s3 .<= ul]
		J4 = J4[ll .<= s4 .<= ul]
	
		s1 = s1[ll .<= s1 .<= ul]
		s2 = s2[ll .<= s2 .<= ul]
		s3 = s3[ll .<= s3 .<= ul]
		s4 = s4[ll .<= s4 .<= ul]




		ax_in = fig.add_axes([0.3, 0.5, 0.3, 0.35])
		ax_in.plot(s1, J1, "r.",ms=6.0)
		ax_in.plot(s2, J2, "r.",ms=6.0)
		ax_in.plot(s3, J3, "r.",ms=6.0)
		ax_in.plot(s4, J4, "r.",ms=6.0)
		ax_in.plot(beta[1:5], Eac[1:5], "bP", ms=8.0)
		#ax_in.plot(beta[1], Eac[1], "bP", ms=8.0)
		#ax_in.plot(beta1[1], Eac1[1], "kP", ms=8.0)
		#ax_in.plot(6.5, 23.1, "bP", ms=8.0)

	
		ax_in.grid(true)
		ax_in.xaxis.set_tick_params(labelsize=36)
		ax_in.yaxis.set_tick_params(labelsize=36)
		@show beta
end
