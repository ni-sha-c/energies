include("../examples/rijke.jl")
include("../src/lss.jl")
function Jac(u,s,n)
	u1 = Rijke_ODE(u,s,n)
	return sum(x->x*x, u1[1:2*Ng])
end
function Jray(u,s,n)
	u1 = Rijke_ODE(u,s,n)
	return dot(cjpixf, u1[1:Ng])
end
function optimize()
nNewton = 5
beta_path = zeros(nNewton)
dJds_mean = zeros(nNewton)

n_samples = 5
n = 2000
dJds = zeros(2,n_samples)
vsh = zeros(N, n, n_samples)
gamma = 0.1
beta_path[1] = 6.5

for nN = 1:nNewton-1
	s = [beta_path[nN], 0.2]
	d = N
	u_trj = zeros(d,n)
	Jac_trj = ones(n)
	Jray_trj = ones(n)
	dJac_trj = zeros(d,n)
	dJray_trj = zeros(d,n)
	f_trj = zeros(d,n)
	X_trj = zeros(d,n)
	nRunup = 60000
	u = Rijke_ODE(rand(d),s,nRunup)
	u_trj[:,end] = u
	for k = 1:n_samples 
		println("sample number:", k)
		u_trj[:,1] = u_trj[:,end]
		X_trj[:,1] = perturbation(u, s)
		# this last parameter is t to be compatible with 
		# ODE problem

		f!(view(f_trj, :, 1), u, s, 1.)

		du_trj = zeros(d, d, n)
		
		du_trj[:,:, 1] = dRijke(u_trj[:,1], s, 1.e-5)
		for i = 2:n
			u_trj[:,i] = Rijke_ODE(u_trj[:,i-1], s, 1)
			vel_p_i = view(u_trj, 1:2*Ng, i)
			Jac_trj[i] = 0.25*sum(x->x*x,vel_p_i)
			Jray_trj[i] = dot(cjpixf, vel_p_i[1:Ng]) 
			dJac_trj[1:2*Ng,i] = 0.5*vel_p_i
			dJray_trj[1:Ng,i] = cjpixf
			X_trj[:,i] = perturbation(u_trj[:,i], s, 1.e-5)
			f!(view(f_trj,:,i), u_trj[:,i], s, 1.)
			du_trj[:,:, i] = dRijke(u_trj[:,i], s, 1.e-5)
		end
		println("set up complete")
		J = [Jac_trj Jray_trj]'
		dJ = reshape([dJac_trj dJray_trj], d, n, 2)
		dJ = permutedims(dJ, [3, 1, 2])
		y, xi = lss(du_trj, X_trj, zeros(d,n), s, 3)
		dJds[:,k] .= compute_sens(y, xi, dJ, f_trj)
		vsh[:,:,k] .= y
	end
	
	thres = 20.0
	# objective function is Jac + alpha*Jray
	alpha = 0.
	n_count = 0
	for i = 1:n_samples
		flag = 0

		for j = 1:n
			if norm(vsh[:,j,i]) > thres
				flag = 1
			end
		end
		if flag == 0
			@show dJds[1,i]
			n_count = n_count + 1
			dJds_mean[nN] += dJds[1,i] + alpha*dJds[2,i]
		end
	end
	dJds_mean[nN] /= n_count
	beta_path[nN+1] = beta_path[nN] - gamma*dJds_mean[nN]
	@show dJds_mean[nN]
end	
return dJds_mean, beta_path, dJds, vsh
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
