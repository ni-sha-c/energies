include("rijke.jl")
	u = zeros(N)
	u[1] = 1.
	s = [7.0, 0.2]
	n = 200
	d = N
	u_trj = zeros(d,n)
	Jac_trj = ones(n)
	Jray_trj = ones(n)
	dJac_trj = zeros(d,n)
	dJray_trj = zeros(d,n)
	xf = 0.2
	cjpixf = cos.(pi*xf.*(1:Ng))
	u_trj[:,1] = u
	for i = 2:n
		u_trj[:,i] = Rijke(u_trj[:,i-1], s, 1)
		vel_p_i = view(u_trj, 1:2*Ng, i)
		Jac_trj[:,i] = 0.25*sum(x->x*x,vel_p_i)
		Jray_trj[:,i] = dot(cjpixf, vel_p_i[1:Ng]) 
		dJac_trj[1:2*Ng,i] = 0.5*vel_p_i
		dJray_trj[1:Ng,i] = cjpixf
	end
	#=
	d = N
	d_u = 3
	dJ = reshape(kron(ones(n),
		[cjpixf; zeros(Ng); zeros(Nc)]),
		d, n)
	println(size(u_trj))
	du_trj = dRijke(u_trj, s, 1.e-5)
	dir_n = zeros(d,n)
	[f!(dir_n[:,i], u_trj[:,i], s, 1.) for i = 1:n]
	X = perturbation(u_trj', s)
	=#
