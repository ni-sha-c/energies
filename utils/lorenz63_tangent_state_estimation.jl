include("../examples/lorenz63.jl")
include("../src/lss.jl")
using LinearAlgebra
using JLD
function refine_parameter(n, n_gd_steps)
	#create epsilon orbit.
	n_runup = 5000
    s = [10., 28., 8/3]
    u0 = rand(3,1)
	u0 = lorenz63(u0, s, n_runup)[end,:,:]
	u_trj = lorenz63(u0, s, n)[:,:,1]
	n, d = size(u_trj)
	eps = 0.0025
	n_noise = 2
	noise = eps*randn(d)
	z_obs = u_trj[:,3]
	u_trj[1,:] .+= noise
	u1 = lorenz63(reshape(u_trj[1,:], d, 1), s, n-1)[:,:,1]
	u_trj = u1
	# run gradient descent to minimize observation error.
	n_gd_steps = 5000
	gamma = 1.0
	dJds = zeros(n_gd_steps)
	d_u = 2
	u = zeros(n, d, n_gd_steps)
	# compute_sens assumes many objective functions
	dJ = zeros(1,d,n)
	error = zeros(n_gd_steps)
	for i = 1:n_gd_steps
		# Set up LSS
		# J is now sum of squared observation error.
		dJ[1,3,:] = -2/n*(z_obs .- u_trj[:,3])
    	du_trj = dlorenz63(u_trj, s)
        X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    	#f = vectorField(u_trj,s)	
    	f = zeros(d,n)
    	y, xi = lss(du_trj, X, f, s, d_u)
		dJds[i] = compute_sens(y, xi, dJ, f)[1]
    	#println(dJds[i])
		# Make u_trj the computed shadowing trajectory
		ds = gamma*dJds[i]
		s[2] = s[2] - ds
		u_trj = u_trj - ds*y'
		error[i] =  sum((u_trj[:,3] .- z_obs).^2)/n
		@show error[i]
		u[:,:,i] = u_trj
	end
	return z_obs, u[:,3,:], error
end
function refine_parameter_and_trajectory(n, n_gd_steps)
	#create epsilon orbit.
	n_runup = 5000
    s = [10., 28., 8/3]
    u0 = rand(3,1)
	u0 = lorenz63(u0, s, n_runup)[end,:,:]
	u_trj = lorenz63(u0, s, n)[:,:,1]
	n, d = size(u_trj)
	eps = 0.001
	n_noise = 2
	noise = eps*randn(d)
	z_obs = u_trj[:,3]
	u_trj[1,:] .+= noise
	u1 = lorenz63(reshape(u_trj[1,:], d, 1), s, n-1)[:,:,1]
	u_trj = u1
	# run gradient descent to minimize observation error.
	#n_gd_steps = 5000
	gamma = 0.1
	dJds = zeros(n_gd_steps)
	d_u = 2
	u = zeros(n, d, n_gd_steps)
	# compute_sens assumes many objective functions
	dJ = zeros(1,d,n)
	error = zeros(n_gd_steps)
	n_pert = 100
	flag = 1
	for i = 1:n_gd_steps
		if flag == 0
			break
		end
		# Set up LSS
		# J is now sum of squared observation error.
		dJ[1,3,:] = -2/n*(z_obs .- u_trj[:,3])
    	du_trj = dlorenz63(u_trj, s)
        X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    	#f = vectorField(u_trj,s)	
    	f = zeros(d,n)
    	y, xi = lss(du_trj, X, f, s, d_u)
		dJds[i] = compute_sens(y, xi, dJ, f)[1]
    	#println(dJds[i])
		# Make u_trj the computed shadowing trajectory
		ds = gamma*dJds[i]
		s[2] = s[2] - ds
		#u_trj = u_trj - ds*y'
		u2 = u_trj[n_pert,:] .- ds*y[:,n_pert]
		u_trj[n_pert:end,:] = lorenz63(reshape(u2,d,1), s, n-n_pert)[:,:,1]
		error[i] =  sum((u_trj[n_pert:end,3] .- z_obs[n_pert:end]).^2)/(n-n_pert)
		@show error[i]
		if i==2
			if error[i] > error[i-1]
				flag = 0
				break
			end
		end
		u[:,:,i] = u_trj
	end
	return z_obs, u[:,3,:], error, flag
end

function assimilate()
	n_repeat = 1
	n = 2000
	n_gd_steps = 5000
	z_obs = zeros(n, n_repeat)
	z_trj = zeros(n, n_gd_steps, n_repeat)
	errors = zeros(n_gd_steps, n_repeat)
	for i = 1:n_repeat
		res1, res2, res3 = refine_parameter(n-1,
											n_gd_steps)    
		z_obs[:,i] = res1
		z_trj[:,:,i] = res2
		errors[:,i] = res3
	end
	println("Errors are ")
	@show errors[:,1]
	save("../data/l63_asmln_thru_param_estn.jld", 
		 "z_obs", z_obs, "z_prd", z_trj, "msq_err",
		 errors)

	return z_obs, z_trj, errors
end
function assimilate_parameter_and_trajectory()
	n_repeat = 100
	n = 2000
	n_gd_steps = 100
	z_obs = zeros(n, n_repeat)
	z_trj = zeros(n, n_gd_steps, n_repeat)
	errors = zeros(n_gd_steps, n_repeat)
	i = 1
	while i <= n_repeat
		@show i
		res1, res2, res3, flag = refine_parameter_and_trajectory(n-1,
											n_gd_steps)    
		if flag==1
			z_obs[:,i] = res1
			z_trj[:,:,i] = res2
			errors[:,i] = res3
			i = i + 1
		end
	end
	println("Errors are ")
	@show errors[:,1]
	save("../data/l63_asmln_ngd200_n2000.jld", 
		 "z_obs", z_obs, "z_prd", z_trj, "msq_err",
		 errors)

	return z_obs, z_trj, errors
end

	
