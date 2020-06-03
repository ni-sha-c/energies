include("solenoid.jl")
include("clvs.jl")
include("lss.jl")
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
	s = [1., 4.]
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
function test_perturbation()
	s = [1.,4.]
	n = 10
	eps = 1.e-4
	u = rand(n,3)
	dFds_fd = (solenoid(u, s .+ eps.*[1,0], 1)[2,:,:] - 
			   solenoid(u, s .- eps.*[1,0], 1)[2,:,:])/
			   (2*eps)
	dFds = perturbation(u, s)
	@test all(isapprox.(dFds_fd - dFds, 0., atol=1.e-8)) == 
	1
end
function test_lss()
	s = [3.,4.]
    m = 1
    n = 1000
    n_runup = 1000
    d = 3
    d_u = 1

	u0 = rand(m,d)
	u_init = transpose(solenoid(u0, s, n_runup
								)[end,:,:])
   	n_samples = 10
	dJds = zeros(n_samples)
	vsh = zeros(d,n+1)
	f = zeros(d,n+1)	
	u = zeros(n+1,d)	
	for i=1:n_samples
		println("Starting LSS, sample ", i)
		u_trj = solenoid(u_init, s, n)[:,:,1]
		du_trj = dsolenoid(u_trj, s)
		X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
		J = objective(u_trj,s)
		dJ = dobjective(u_trj, s) 
		v, dJds[i] = lss(u_trj, du_trj, X, f, J, dJ, s, d_u)
		vsh .= v
		u .= u_trj
		u_init .= reshape(u_trj[end,:],1,3)
		println(dJds[i])	
	end
	# check with analytical shadowing direction
	# which is [r,0,0]
	n_spinup = 100
	u = u[n_spinup+1:end,:]
	vsh = vsh[:, n_spinup+1:end]
	n = n - n_spinup
	x, y, z = view(u,:,1), view(u,:,2),
				  view(u,:,3)
	r = sqrt.(x.*x .+ y.*y .+ z.*z)
	t_dir = [-y./r  x./r  zeros(n+1)]'
	r_dir = [x./r   y./r  zeros(n+1)]'
	prj_theta = [dot(t_dir[:,i], vsh[:,i]) for i=1:n] 
	prj_z = [dot([0,0,1.], vsh[:,i]) for i=1:n] 
	#@test all(isapprox.(prj_theta, 0., atol=0.01)) == 1
	#@test all(isapprox.(prj_z, 0., atol=0.01)) == 1
	@test isapprox((sum(dJds)/n_samples),1.0,rtol=0.1)  
	return vsh, dJds
end
