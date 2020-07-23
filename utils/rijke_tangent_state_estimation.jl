include("../examples/rijke.jl")
include("../src/lss.jl")
using LinearAlgebra
using JLD
using OrdinaryDiffEq
function refine_parameter_and_trajectory(n, n_gd_steps)
    #create epsilon orbit.
    n_runup = 400000
    s = [7.0,0.2]
    d = N
    u0 = Rijke_ODE(rand(d), s, n_runup)
    eps = 0.05
    noise = eps*randn(d)
    u_trj = zeros(d,n)
    u_trj[:,1] = u0 + noise
    u_obs = zeros(d,n)
    u_obs[:,1] = u0
    for i = 2:n
        u_trj[:,i] = Rijke_ODE(u_trj[:,i-1], s, 1)
        u_obs[:,i] = Rijke_ODE(u_obs[:,i-1], s, 1)
    end
    z_prd = zeros(n, n_gd_steps)
    z_obs = sum(u_obs[1:2*Ng,:].*u_obs[1:2*Ng,:],dims=1)./4    
    z_trj = sum(u_trj[1:2*Ng,:].*u_trj[1:2*Ng,:],dims=1)./4
    z_prd[:,1] = z_trj
    # run gradient descent to minimize observation error.
    gamma = 0.1
    dJds = zeros(n_gd_steps)
    d_u = 2
    # compute_sens assumes many objective functions
    dJ = zeros(1,d,n)
    error = zeros(n_gd_steps)
    n_pert = 100
    flag = 1

    du_trj = zeros(d,d,n)
    f_trj = zeros(d,n)
    X_trj = zeros(d,n)
    @time for i = 1:n_gd_steps
        if flag == 0
            break
        end
        # Set up LSS
        # J is now mean of squared observation error.
        for k = 1:n
            z_prd[k,i] = sum(x->x*x, u_trj[1:2*Ng,k])/4
            error_k = (z_obs[k] .- z_prd[k,i])^2.0/n
            dJ[1,1:2*Ng,k] = -1/n*error_k*u_trj[1:2*Ng,k]
            du_trj[:,:,k] = dRijke(u_trj[:,k], s, 1.e-6)
            X_trj[:,k] = perturbation(u_trj[:,k],s,1.e-6) #ith col in T_{u_{i+1}} M
            f!(view(f_trj,:,k), u_trj[:,k], s, 1.)    
            #f = zeros(d,n)
        end
        y, xi = lss(du_trj, X_trj, f_trj, s, d_u)
        dJds[i] = compute_sens(y, xi, dJ, f_trj)[1]
        #println(dJds[i])
        # Make u_trj the computed shadowing trajectory
        ds = gamma*dJds[i]
        s[1] = s[1] - ds
        #u_trj = u_trj - ds*y'
        u2 = u_trj[:,n_pert] .- ds*y[:,n_pert]
        for k = n_pert+1:n
            u_trj[:,k] = Rijke_ODE(u_trj[:,k-1], s, 1)
            error[i] += (z_prd[k,i] - z_obs[k])^2.0/(n-n_pert)
        end
        @show error[i]
        if i==2
            if error[i] > error[i-1]
                flag = 0
                break
            end
        end
    end
    return z_obs, z_prd, error, flag
end
function assimilate_parameter_and_trajectory(filename)
    n_repeat = 1
    n = 2000
    n_gd_steps = 200
    z_obs = zeros(n, n_repeat)
    z_trj = zeros(n, n_gd_steps, n_repeat)
    errors = zeros(n_gd_steps, n_repeat)
    i = 1
    while i <= n_repeat
        @show i
        res1, res2, res3, flag = refine_parameter_and_trajectory(n,
                                            n_gd_steps)    
        if flag==1
            z_obs[:,i] = res1
            z_trj[:,:,i] = res2
            errors[:,i] = res3
            i = i + 1
        end
    end
    #println("Errors are ")
    #@show errors[:,1]
    save(filename, 
         "z_obs", z_obs, "z_prd", z_trj, "msq_err",
         errors)

    #return z_obs, z_trj, errors
end

    
