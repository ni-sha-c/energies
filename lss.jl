include("lorenz63.jl")
include("clvs.jl")
"""
clvs : nxdxd_u
"""
	s = [10., 28., 8/3]
	m = 1
	n = 25000
	n_runup = 5000
	u0 = rand(3,m)
	u_init = lorenz63(u0, s, n)[end,:,:]
	u_trj = lorenz63(u_init, s, n)[:,:,1]
	du_trj = dlorenz63(u_trj, s)
	du_trj = permutedims(du_trj,[2,3,1])
	les, clv_trj = clvs(du_trj,1)

