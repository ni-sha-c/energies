include("lorenz63.jl")
include("clvs.jl")
using Test
function test_dlorenz63()
	s = [10., 28., 8/3]
	m = 10
	u0 = rand(3,m)
	u_trj = lorenz63(u0, s, 1)[2,:,:]
	eps = 1.e-4
	du_x = (lorenz63(u_trj .+ eps.*[1; 0.; 0.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[1; 0.; 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_y = (lorenz63(u_trj .+ eps.*[0; 1.; 0.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[0; 1.; 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_z = (lorenz63(u_trj .+ eps.*[0; 0.; 1.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[0; 0.; 1.], 
						s, 1)[2,:,:])./(2*eps)
	du_fd = reshape(collect([du_x' du_y' du_z']), 
						m, 3, 3)
	du = dlorenz63(u_trj', s)
	@assert all(isapprox.(du_fd, du)) == 1
end
function test_les()
	s = [10., 28., 8/3]
	m = 1
	n = 25000
	u0 = rand(3,m)
	u_trj = lorenz63(u0, s, n)[:,:,1]
	du_trj = dlorenz63(u_trj, s)
	du_trj = permutedims(du_trj,[2,3,1])
	les, clv_trj = clvs(du_trj,3)
	println(les/0.005)
	@assert all(isapprox.(les/0.005, [0.91, 0., 
		-14.572], rtol=1.e-1)) == 1
end
function test_lss()
	s = [10., 28., 8/3]
    m = 1
    n = 2000
    n_runup = 5000
    u0 = rand(3,m)
    u_init = lorenz63(u0, s, n)[end,:,:]
    d = 3
    d_u = 1
	dJ = zeros(d,n+1)
	dJ[3,:] .= 1.
	include("lss.jl")
	n_samples = 10
	dJds = zeros(n_samples)
	for i=1:n_samples
		println("Starting LSS, sample ", i)
		u_trj = lorenz63(u_init, s, n)[:,:,1]
		du_trj = dlorenz63(u_trj, s)
    	du_trj = permutedims(du_trj,[2,3,1])
    	X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
		f = vectorField(u_trj,s)	
		J = u_trj[:,3]
		vsh, dJds[i] = lss(u_trj, du_trj, X, f, J, dJ, 
						  s, d_u)

		global u_init .= reshape(u_trj[end,:],3,1)
	end
	@test isapprox((sum(dJds)/n_samples),1.0,rtol=0.1)  
	return vsh, dJds
end
