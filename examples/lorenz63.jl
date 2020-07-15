dt = 0.005
function lorenz63(u0, s, n)
    sigma, rho, beta = s
    d, m = size(u0)
    n = n+1
    u_trj = zeros((m,d,n))
    u_trj[:,:,1] = u0'
    for i = 2:n
    	x = u_trj[:,1,i-1]
    	y = u_trj[:,2,i-1]
    	z = u_trj[:,3,i-1]
    
    	u_trj[:,1,i] = x + dt*(sigma.*(y - x))
    	u_trj[:,2,i] = y + dt*(x.*(rho .- z) - y)
    	u_trj[:,3,i] = z + dt*(x.*y - beta.*z)
    end 
    return permutedims(u_trj,[3,2,1])
end
function dlorenz63(u, s)
    sigma, rho, beta = s
    n, d = size(u)
    x = view(u,:,1)
    y = view(u,:,2)
    z = view(u,:,3)
    du = zeros(n, d, d)
    @. du[:,1,1] = 1.0 - dt*sigma
    @. du[:,1,2] = dt*sigma
    @. du[:,2,1] = dt*(rho - z) 
    @. du[:,2,2] = 1.0 - dt
    @. du[:,2,3] = -dt*x 
    @. du[:,3,1] = dt*y
    @. du[:,3,2] = dt*x
    @. du[:,3,3] = 1.0 - dt*beta
    return reshape([du[:,:,1]'; du[:,:,2]'; 
    				du[:,:,3]'], d, d, n)
end
function perturbation(u,s)
    n, d = size(u)
    # the perturbation in row i in T_{u_(i+1)} M
    return [zeros(1,n); dt*u[:,1]'; zeros(1,n)]
end
function vectorField(u,s)
    n, d = size(u)
    sigma, rho, beta = s
    u = u'
    x, y, z = u[1,:], u[2,:], u[3,:]
    return [sigma.*(y - x)  x.*(rho .- z) - y  x.*y - beta.*z]'
end
function lorenz63_rhs_ad(du, u, s, t)
    du[1] = s[1]*(u[2] - u[1])
	du[2] = u[1]*(s[2] - u[3]) - u[2]
	du[3] = u[1]*u[2] - s[3]*u[3]
end
function lorenz63_ad(u0, s, n)
	t = n*dt
    prob = ODEProblem(lorenz63_rhs_ad, u0, (0.,t), s)
	sol = Array(solve(prob, Tsit5(), saveat=dt))
	return sol[:,end]
end
function obj_fun(u0, s)
    prob = ODEProblem(lorenz63_rhs_ad, u0, (0.,1.1), s)
	#_prob = remake(prob,u0=u0,p=s) 
	sol = solve(prob, Tsit5(), saveat=0.005)
	sum(sol[3,:])/size(sol,2)
end
