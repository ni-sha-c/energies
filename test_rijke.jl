include("rijke.jl")
include("lss.jl")
	s = [7.0, 0.2]
	n = 5000
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
	dJds = zeros(2,n_samples)
	vsh = zeros(d, n, n_samples)
	nRunup = 40000
	u = Rijke(rand(d),s,nRunup)
	u_trj[:,end] = u
	for k = 1:n_samples 
		println("sample number:", k)
		u_trj[:,1] = u_trj[:,end]
		X_trj[:,1] = perturbation(u, s)
		# this last parameter is t to be compatible with 
		# ODE problem

		f!(view(f_trj, :, 1), u, s, 1.)

		du_trj = zeros(d, d, n)
		
		du_trj[:,:, 1] = dRijke(u_trj[:,1], s, 1.e-6)
		for i = 2:n
			u_trj[:,i] = Rijke(u_trj[:,i-1], s, 1)
			vel_p_i = view(u_trj, 1:2*Ng, i)
			Jac_trj[i] = 0.25*sum(x->x*x,vel_p_i)
			Jray_trj[i] = dot(cjpixf, vel_p_i[1:Ng]) 
			dJac_trj[1:2*Ng,i] = 0.5*vel_p_i
			dJray_trj[1:Ng,i] = cjpixf
			X_trj[:,i] = perturbation(u_trj[:,i], s, 1.e-6)
			f!(view(f_trj,:,i), u_trj[:,i], s, 1.)
			du_trj[:,:, i] = dRijke(u_trj[:,i], s, 1.e-6)
		end
		J = [Jac_trj Jray_trj]'
		dJ = reshape([dJac_trj dJray_trj], d, n, 2)
		dJ = permutedims(dJ, [3, 1, 2])
		y, xi = lss(du_trj, X_trj, f_trj, s, 3)
		dJds[:,k] = compute_sens(y, xi, dJ, f_trj)
		vsh[:,:,k] = y
	end
