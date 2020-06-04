include("rijke.jl")
	u = zeros(N)
	u[1] = 1.
	s = [7.0, 0.2]
	n = 200
	d = N
	u_trj = zeros(d,n)
	J_trj = -ones(n)
	u_trj[:,1] = u
	for i = 2:n
		u_trj[:,i] = Rijke(u_trj[:,i-1], s, 1) 
		J_trj[i] = J[1]
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
