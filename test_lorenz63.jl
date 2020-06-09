include("lorenz63.jl")
include("clvs.jl")
include("lss.jl")
using Test
#using DiffEqSensitivity, OrdinaryDiffEq, Zygote
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
    du_fd = reshape(collect([du_x; du_y; du_z]), 
    					3, 3, m)
    du = dlorenz63(u_trj', s)
    @test all(isapprox.(du_fd - du, 0.,atol=1.e-8)) == 1
end
function test_les()
    s = [10., 28., 8/3]
    m = 1
    n = 20000
    u0 = rand(3,m)
    u_trj = lorenz63(u0, s, n)[:,:,1]
    du_trj = dlorenz63(u_trj, s)
    les, clv_trj = clvs(du_trj,3)
    println(les/0.005)
    @test all(isapprox.(les/0.005, [0.91, 0., 
    	-14.572], rtol=1.e-1)) == 1
end
function test_lss()
    s = [10., 28., 8/3]
    m = 1
    n = 5000
    n_runup = 5000
    u0 = rand(3,m)
    u_init = lorenz63(u0, s, n)[end,:,:]
    d = 3
    d_u = 2
    dJ = zeros(d,n+1)
    dJ[3,:] .= 1.
    n_samples = 100
    dJds = zeros(n_samples)
    vsh = zeros(d,n+1,n_samples)
    for i=1:n_samples
    	println("Starting LSS, sample ", i)
    	u_trj = lorenz63(u_init, s, n)[:,:,1]
    	du_trj = dlorenz63(u_trj, s)
        X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    	#f = vectorField(u_trj,s)	
    	f = zeros(d,n+1)
    	J = u_trj[:,3]
    	y, dJds[i] = lss(u_trj, du_trj, X, f, J, 
    					 dJ, s, d_u)
    	println(dJds[i])
    	vsh[:,:,i] = y
    	u_init .= reshape(u_trj[end,:],3,1)
    end
    @test isapprox((sum(dJds)/n_samples),1.0,rtol=0.1)  
    return vsh, dJds
end
function test_condition_number()
    s = [10., 28., 8/3]
    m = 1
    n_runup = 5000
    u0 = rand(3,m)
    u_init = lorenz63(u0, s, n_runup)[end,:,:]
    d = 3
    d_u = 2
    n_samples = 30
    dJds = zeros(n_samples)
    condnum = zeros(n_samples)
    n_arr = LinRange(500, 5000, n_samples)
    for i=1:n_samples
    	println("Starting LSS, sample ", i)
    	n = floor(Int64,n_arr[i])
    	u_trj = lorenz63(u_init, s, n)[:,:,1]
    	du_trj = dlorenz63(u_trj, s)
        X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    	f = zeros(d,n+1)
    	J = u_trj[:,3]
    	dJ = zeros(d,n+1)
    	dJ[3,:] .= 1.

    	y, dJds[i], c = lss(u_trj, du_trj, X, f, J, 
    					 dJ, s, d_u)
    	condnum[i] = c
    	println(dJds[i])
    	u_init .= reshape(u_trj[end,:],3,1)
    end
    @test isapprox((sum(dJds)/n_samples),1.0,rtol=0.1)  
    return dJds, condnum
end
function test_Zygote()

    s = [10., 28., 8/3]
	du0, ds = zeros(3), 0.
	n_samples = 10000
	#prob = ODEProblem(lorenz63_ad, u0, (0.,1.), s) 
    for i =1:n_samples
		u0 = rand(3)
		println("sample number, ", i)
		#sol = solve(prob,Tsit5(),saveat=0.005,sensealg=QuadratureAdjoint())
    	du01, ds1 = Zygote.gradient(obj_fun, u0, s)	
		du0 += du01/n_samples
		ds += ds1/n_samples
	end
end

