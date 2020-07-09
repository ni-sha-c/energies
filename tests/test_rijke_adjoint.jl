include("../examples/rijke.jl")
include("../src/adjoint_lss.jl")
function Jac(u,s,n)
	u1 = Rijke_ODE(u,s,n)[:,end]
	return sum(x->x*x, u1[1:2*Ng])
end
function Jray(u,s,n)
	u1 = Rijke_ODE(u,s,n)[:,end]
	return dot(cjpixf, u1[1:Ng])
end
	s = [6.5, 0.2]
	n = 1000
	d = N
	u_trj = zeros(d,n)
	Jac_trj = ones(n)
	Jray_trj = ones(n)
	dJac_trj = zeros(d,n)
	dJray_trj = zeros(d,n)
	f_trj = zeros(d,n)
	X_trj = zeros(d,n)
	n_samples = 1
	dJds = zeros(2,n_samples)
	vsh = zeros(d, n, n_samples)
	nRunup = 1000000
	u = Rijke_ODE(rand(d),s,nRunup)
	u_trj[:,end] = u
	#for k = 1:n_samples 
		println("sample number:", k)
		u_trj[:,1] = u_trj[:,end]
		X_trj[:,1] = perturbation(u, s)
		# this last parameter is t to be compatible with 
		# ODE problem

		f!(view(f_trj, :, 1), u, s, 1.)

		du_trj = zeros(d, d, n)
		
		du_trj[:,:, 1] = dRijke(u_trj[:,1], s, 1.e-6)
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
		println("Reversing time...")
		
		dJ = reverse(dJac_trj, dims=2)
		du_trj = reverse(du_trj, dims=3)
		du_trj = permutedims(du_trj, [2, 1, 3])
		X_trj = reverse(X_trj, dims=2)
		X_trj = reshape(X_trj, 1, d, n)
		f_trj = reverse(f_trj, dims=2)
		u_next = Rijke_ODE(u_trj[:,end], s, 1)
		f_end = zeros(d)
		f!(f_end, u_next, s, 1.)
		f_trj = [f_end f_trj[:,1:end-1]]
		y, xi = lss(du_trj, dJ, zeros(d,n), s, 1, f_trj)
		dJds[1,k] = compute_sens(y, xi, X_trj, zeros(d,n))[1]
		vsh[:,:,k] = y
	#end
	println("Shadowing complete. Ensemble sensitivity starting...")
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
