include("../examples/lorenz63.jl")
include("../src/lss.jl")
using LinearAlgebra
function refine_parameter(n)
	#create epsilon orbit.
	n_runup = 5000
    s = [10., 28., 8/3]
    u0 = rand(3,1)
	u0 = lorenz63(u0, s, n_runup)[end,:,:]
	u_trj = lorenz63(u0, s, n)[:,:,1]
	n, d = size(u_trj)
	eps = 0.1
	noise = eps*randn(n)
	u_obs = copy(u_trj)
	u_obs[:,3] .+= noise
	# run gradient descent to minimize observation error.
	n_gd_steps = 4
	gamma = 0.1
	dJds = zeros(n_gd_steps)
	d_u = 2
	u = zeros(n, d, n_gd_steps)
	# compute_sens assumes many objective functions
	dJ = zeros(1,d,n)
	error = zeros(n_gd_steps)
	for i = 1:n_gd_steps
		# Set up LSS
		# J is now sum of squared observation error.
		dJ[1,3,:] = -2.0*noise
    	du_trj = dlorenz63(u_trj, s)
        X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    	#f = vectorField(u_trj,s)	
    	f = zeros(d,n)
    	y, xi = lss(du_trj, X, f, s, d_u)
		dJds[i] = compute_sens(y, xi, dJ, f)[1]
    	println(dJds[i])
		# Make u_trj the computed shadowing trajectory
		ds = gamma*dJds[i]
		s[2] = s[2] - ds
		u_trj = u_trj + ds*y'
		error[i] =  sum((u_trj[:,3] .- u_obs[:,3]).^2)/n 
		u[:,:,i] = u_trj
	end
	return u_obs[:,3], u[:,3,:], error
end
function assimilate()
	n_repeat = 1000
	n = 2000
	n_gd_steps = 4
	z_obs = zeros(n, n_repeat)
	z_trj = zeros(n, n_gd_steps, n_repeat)
	errors = zeros(n_gd_steps, n_repeat)
	for i = 1:n_repeat
		res1, res2, res3 = refine_parameter(n-1)    
		z_obs[:,i] = res1
		z_trj[:,:,i] = res2
		errors[:,i] = res3
	end
	return z_obs, z_trj, errors
end

	
