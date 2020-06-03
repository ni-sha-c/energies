include("solenoid.jl")
include("clvs.jl")
using Test
function test_dsolenoid()
	s = [1.,4]
	m = 10
	u_trj = rand(m,3)
	eps = 1.e-5
	du_x = (solenoid(u_trj .+ eps.*[1 0. 0.],
						s, 1)[2,:,:] - 
				solenoid(u_trj .- eps.*[1 0. 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_y = (solenoid(u_trj .+ eps.*[0 1. 0.],
						s, 1)[2,:,:] - 
				solenoid(u_trj .- eps.*[0 1. 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_z = (solenoid(u_trj .+ eps.*[0 0. 1.],
						s, 1)[2,:,:] - 
				solenoid(u_trj .- eps.*[0 0. 1.], 
						s, 1)[2,:,:])./(2*eps)
	du_fd = reshape(collect([du_x; du_y; du_z]), 
						3, 3, m)
	du = dsolenoid(u_trj, s)
	@assert all(isapprox.(du_fd, du, rtol=1.e-2)) == 1
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
	n_samples = 1
	dJds = zeros(n_samples)
	vsh = zeros(d,n+1)
	for i=1:n_samples
		println("Starting LSS, sample ", i)
		u_trj = lorenz63(u_init, s, n)[:,:,1]
		du_trj = dlorenz63(u_trj, s)
    	du_trj = permutedims(du_trj,[2,3,1])
    	X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
		f = vectorField(u_trj,s)	
		J = u_trj[:,3]
		y, dJds[i] = lss(u_trj,  
						du_trj, X, f, J, dJ, 
						  s, d_u)
		println(size(y))
		vsh .= y
		u_init .= reshape(u_trj[end,:],3,1)
	end
	@test isapprox((sum(dJds)/n_samples),1.0,rtol=0.1)  
	return vsh, dJds
end
