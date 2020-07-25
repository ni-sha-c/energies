include("../examples/rijke.jl")
include("../src/lss.jl")
using JLD
using OrdinaryDiffEq
function Rijke_tangent_sensitivity(n_spe)
	s = [7.0, 0.2]
	n = 2000
	d = N

	du_trj = zeros(d, d, n)
	dJ_trj = zeros(2,d,n)
	f_trj = zeros(d,n)
	X_trj = zeros(d,n)
	
	dJds = zeros(2)
	vsh = zeros(d, n)
	nRunup = 1000000
	u = Rijke_ODE(rand(d),s,nRunup)

	for i = 1:n
		vel_p_i = view(u, 1:2*Ng)
		dJ_trj[1,1:2*Ng,i] = 0.5*vel_p_i
		dJ_trj[2,1:Ng,i] = cjpixf
		X_trj[:,i] = perturbation(u, s, 1.e-6)
		f!(view(f_trj,:,i), u, s, 1.)
		du_trj[:,:, i] = dRijke(u, s, 1.e-6)
		u = Rijke_ODE(u, s, 1)
	end
	y, xi = lss(du_trj, X_trj, f_trj, s, 2)
	dJds = compute_sens(y, xi, dJ_trj, f_trj)
	vsh[:,:] = y
	save(string("../data/rijke_tangent_sensitivity/",
		    "vsh_and_dJds_", string(n_spe), 
		    ".jld"), "vsh", vsh, "dJds", dJds)
end
function collect_sensitivities()
	pmap(Rijke_tangent_sensitivity, 241:480)
end

#=
	n_samples = 1000
	dJacds_ens, dJrayds_ens = zeros(2,n_samples), 
							zeros(2,n_samples)
	up = zeros(d)
	um = zeros(d)
	eps = 1.e-3
	n = 100 # produces a growth of ≈ e
	for i = 1:n_samples
			@time dJacds_ens[:,i] = Zygote.gradient(s ->Jac(u, s, n), s)[1]
		dJrayds_ens[:,i] = Zygote.gradient(s ->Jray(u, s, n),
										   s)[1]
		u .= Rijke(u, [7.0, 0.2], n)
	end
=#
