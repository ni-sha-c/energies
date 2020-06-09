include("rijke.jl")
	s = [7.0, 0.2]
	n = 200
	d = N
	u_trj = zeros(d,n)
	Jac_trj = ones(n)
	Jray_trj = ones(n)
	dJac_trj = zeros(d,n)
	dJray_trj = zeros(d,n)
	f_trj = zeros(d,n)
	X_trj = zeros(d,n)
	xf = 0.2
	cjpixf = cos.(pi*xf.*(1:Ng))
	n_samples = 1
	dJac_ds = zeros(n_samples)
	vsh = zeros(d, n_samples)
	for k = 1:n_samples 
		u = zeros(N)
		u[1] = 1.

		u_trj[:,1] = u
		X_trj[:,1] = perturbation(u, s)
		# this last parameter is t to be compatible with 
		# ODE problem

		f!(view(f_trj, :, 1), u, s, 1.) 
		du_trj = zeros(d, d, n)
		for i = 2:n
			u_trj[:,i] = Rijke(u_trj[:,i-1], s, 1)
			vel_p_i = view(u_trj, 1:2*Ng, i)
			Jac_trj[i] = 0.25*sum(x->x*x,vel_p_i)
			Jray_trj[i] = dot(cjpixf, vel_p_i[1:Ng]) 
			dJac_trj[1:2*Ng,i] = 0.5*vel_p_i
			dJray_trj[1:Ng,i] = cjpixf
			X_trj[:,i] = perturbation(u_trj[:,i], s)
			f!(view(f_trj,:,i), u_trj[:,i], s, 1.)
			du_trj[:, i] = dRijke(u_trj[:,i], s, 1.e-5)
		end
		y, dJac_ds[k] = lss(u_trj, du_trj, X_trj, 
							f_trj, Jac_trj, dJac_trj,
							s, 3)
		vsh[:,k] = y
	end
