include("../examples/rijke.jl")
include("../src/clvs.jl")
#using PyPlot
using OrdinaryDiffEq
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function compute_acoustic_energy()
    n_samples = 112
    d_u = 3
    Eac = SharedArray{Float64}(n_samples)
    les = SharedArray{Float64}(d_u, n_samples)
    les .= 0.
    Eac .= 0.
    beta = LinRange(2.0,9.0,n_samples)
    nSteps = 5000
    nRunup = 400000
    answer = @distributed for i = 1:n_samples
    println("beta = ", beta[i])
    s = [beta[i], 0.2]
    u = Rijke_ODE(rand(N), s, nRunup)
    du_trj = zeros(N, N, nSteps)
    for j = 1:nSteps    
        u = Rijke_ODE(u, s, 1)
        du_trj[:,:,j] = dRijke(u, s, 1.e-6)
        Eac[i] += norm(u[1:2*Ng])^2.0/4/nSteps
    end
    les[:,i], clv_trj = clvs(du_trj, d_u) 
    @show Eac[i]
    end
    wait(answer)
    save("../data/attractor/Eac1.jld", "beta", beta, 
     "Eac", Eac, "les", les)
    return beta, Eac, les

end
function plot_acoustic_energy_colored()
	X = load("../data/attractor/more_les/Eac.jld")
	Y = load("../data/attractor/more_les/Eac1.jld")
	Z = load("../data/attractor/Eac.jld")
	W = load("../data/attractor/Eac_and_beta.jld")
	Eac1 = X["Eac"]
	Eac2 = Y["Eac"]
	beta1 = X["beta"]
	beta2 = Y["beta"]
	les1 = X["les"]/dt
	les2 = Y["les"]/dt
	Eac3 = Z["Eac"]
	Eac4 = W["Eac"]
	beta3 = Z["beta"]
	beta4 = W["beta"]

	b_c_ll = 6.2
	b_c_ul = 7.4
	b_qp_ll_1 = 5.0	
	b_qp_ul_1 = 6.2
	b_qp_ll_2 = 7.4
	b_qp_ul_2 = 8.0

	beta1_c = beta1[b_c_ll .< beta1 .< b_c_ul]
	chaotic_Eac1 = Eac1[b_c_ll .< beta1 .< b_c_ul]
	beta2_c = beta2[b_c_ll .< beta2 .< b_c_ul]
	chaotic_Eac2 = Eac2[b_c_ll .< beta2 .< b_c_ul]
	
	beta3_c = beta3[b_c_ll .< beta3 .< b_c_ul]
	chaotic_Eac3 = Eac3[b_c_ll .< beta3 .< b_c_ul]
	beta4_c = beta4[b_c_ll .< beta4 .< b_c_ul]
	chaotic_Eac4 = Eac4[b_c_ll .< beta4 .< b_c_ul]
	
	beta1_qp_1 = beta1[b_qp_ll_1 .< beta1 .< b_qp_ul_1]
	qp1_Eac1 = Eac1[b_qp_ll_1 .< beta1 .< b_qp_ul_1]
	beta2_qp_1 = beta2[b_qp_ll_1 .< beta2 .< b_qp_ul_1]
	qp1_Eac2 = Eac2[b_qp_ll_1 .< beta2 .< b_qp_ul_1]
	
	beta1_qp_2 = beta1[b_qp_ll_2 .< beta1 .< b_qp_ul_2]
	qp2_Eac1 = Eac1[b_qp_ll_2 .< beta1 .< b_qp_ul_2]
	beta2_qp_2 = beta2[b_qp_ll_2 .< beta2 .< b_qp_ul_2]
	qp2_Eac2 = Eac2[b_qp_ll_2 .< beta2 .< b_qp_ul_2]

	beta3_qp_1 = beta3[b_qp_ll_1 .< beta3 .< b_qp_ul_1]
	qp1_Eac3 = Eac3[b_qp_ll_1 .< beta3 .< b_qp_ul_1]
	beta4_qp_1 = beta4[b_qp_ll_1 .< beta4 .< b_qp_ul_1]
	qp1_Eac4 = Eac4[b_qp_ll_1 .< beta4 .< b_qp_ul_1]

	beta3_qp_2 = beta3[b_qp_ll_2 .< beta3 .< b_qp_ul_2]
	qp2_Eac3 = Eac3[b_qp_ll_2 .< beta3 .< b_qp_ul_2]
	beta4_qp_2 = beta4[b_qp_ll_2 .< beta4 .< b_qp_ul_2]
	qp2_Eac4 = Eac4[b_qp_ll_2 .< beta4 .< b_qp_ul_2]


	fig, ax = subplots(1,1)
    ax.plot(beta1, Eac1, "b.")
    ax.plot(beta2, Eac2, "b.")
    ax.plot(beta1_c, chaotic_Eac1, "r.")
    ax.plot(beta2_c, chaotic_Eac2, "r.")
	ax.plot(beta1_qp_1, qp1_Eac1, "g.")
	ax.plot(beta1_qp_2, qp2_Eac1, "g.")
	ax.plot(beta2_qp_1, qp1_Eac2, "g.")
	ax.plot(beta2_qp_2, qp2_Eac2, "g.")
	ax.plot(beta3, Eac3, "b.")
    ax.plot(beta4, Eac4, "b.")
    ax.plot(beta3_c, chaotic_Eac3, "r.")
    ax.plot(beta4_c, chaotic_Eac4, "r.")
	ax.plot(beta3_qp_1, qp1_Eac3, "g.")
	ax.plot(beta3_qp_2, qp2_Eac3, "g.")
	ax.plot(beta4_qp_1, qp1_Eac4, "g.")
	ax.plot(beta4_qp_2, qp2_Eac4, "g.")



    ax.xaxis.set_tick_params(labelsize=25)
    ax.yaxis.set_tick_params(labelsize=25)
    ax.set_xlabel(L"$\beta$",fontsize=25)
    ax.set_ylabel(L"$<E_{\rm ac}>$",fontsize=25)
	grid(true)
	#=
    fig1, ax1 = subplots(1,1)
	ax1.plot(beta1, les1[1,:], ".-")
	ax1.plot(beta1, les1[2,:], ".-")
	ax1.plot(beta1, les1[3,:], ".-")
    ax1.xaxis.set_tick_params(labelsize=25)
    ax1.yaxis.set_tick_params(labelsize=25)
    ax1.set_xlabel(L"$\beta$",fontsize=25)
    ax1.set_ylabel(L"$\lambda$",fontsize=25)
	grid(true)

    fig2, ax2 = subplots(1,1)
	ax2.plot(beta2, les2[1,:], ".-")
	ax2.plot(beta2, les2[2,:], ".-")
	ax2.plot(beta2, les2[3,:], ".-")
    ax2.xaxis.set_tick_params(labelsize=25)
    ax2.yaxis.set_tick_params(labelsize=25)
    ax2.set_xlabel(L"$\beta$",fontsize=25)
    ax2.set_ylabel(L"$\lambda$",fontsize=25)
	grid(true)
	=#
end
