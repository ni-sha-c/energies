include("../examples/rijke.jl")
include("../src/lss.jl")
using OrdinaryDiffEq
using JLD
using SharedArrays
function Rijke_tangent_sensitivity(beta)
	s = [beta, 0.2]
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
	return dJds[1]
end
function average_sensitivity(beta)
	n_samples = 100
	dJds_avg = @distributed (+) for i = 1:n_samples
        Rijke_tangent_sensitivity(beta)
	end
	dJds_avg /= n_samples
end
function optimize()
    nNewton = 20
    beta_path = zeros(nNewton)
    dJds = zeros(nNewton)
    gamma = 0.1
    beta_path[1] = 6.5
    beta = beta_path[1]
    for nN = 1:nNewton
	dJds[nN] = average_sensitivity(beta)
       	beta = beta - gamma*dJds[nN]
	beta_path[nN] = beta
    end
    save("../data/param_optim/beta_dJds_1.jld",
	 "beta_path", beta_path,
	 "dJds_path", dJds)
end
function compute_Eac(beta)
	nRunup = 1000000
	u = Rijke_ODE(rand(d),s,nRunup)
	nSteps = 20000
	Eac = 0.
	for i = 1:nSteps
		u = Rijke_ODE(u, s, 1)
		Eac += 0.25*norm(u[1:2*Ng])^2.0/nSteps
	end
	return Eac
end
function compute_Eac_along_path()
	X = load(string("../data/param_optim/",
			"beta_dJds_1.jld"))
	Y = load(string("../data/param_optim/",
			"beta_dJds.jld"))
	beta1 = X["beta_path"]
	n1, =  size(beta1)
	beta = Y["beta_path"]
	beta = [beta1; beta]
	n, = size(beta)
	Eac = SharedArray{Float64}(n)
	@distributed for i=1:n
		Eac[i] = compute_Eac(beta[i])
	end
	save("../data/param_optim/Eac_along_path.jld",
		"Eac1", Eac[1:n1], 
		"Eac", Eac[n1+1:end],
		"beta", beta,
		"beta1", beta1)
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
