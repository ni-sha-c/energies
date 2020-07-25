include("../examples/rijke.jl")
include("../src/adjoint_lss.jl")
using JLD
using OrdinaryDiffEq
function Rijke_adjoint_sensitivity(n_spe)
    s = [7.0, 0.2]
    n = 2000
    d = N

    J_trj = ones(n)
    dJ_trj = zeros(d,n)
    f_trj = zeros(d,n)
    X_trj = zeros(2,d,n)
    du_trj = zeros(d, d, n)
    
    dJds = zeros(2)
    vsh = zeros(d, n)
    nRunup = 1000000
    u = Rijke_ODE(rand(d),s,nRunup)
    
    for i = 1:n
        du_trj[:,:, i] = dRijke(u, s, 1.e-6)
        X_trj[1,:,i] = perturbation(u, s)
        X_trj[2,:,i] = tau_perturbation(u, s)
           
        #this last parameter is t to be compatible with 
        # ODE problem
        f!(view(f_trj, :, i), u, s, 1.)

        vel_p_i = view(u, 1:2*Ng)
        J_trj[i] = 0.25*sum(x->x*x,vel_p_i)
        dJ_trj[1:2*Ng,i] = 0.5*vel_p_i
    
        u = Rijke_ODE(u, s, 1)

    end
    du_trj = reverse(du_trj, dims=3)
    du_trj = permutedims(du_trj, [2, 1, 3])


    u = Rijke_ODE(u, s, 1)
    X_trj = reverse(X_trj, dims=3)
    X1_end = perturbation(u, s)
    X2_end = tau_perturbation(u, s)
    X_trj[1,:,:] = [X1_end X_trj[1,:,1:end-1]]
    X_trj[2,:,:] = [X2_end X_trj[2,:,1:end-1]]


    dJ_trj = reverse(dJ_trj, dims=2)
    dJ_trj /= n

    f_trj = reverse(f_trj, dims=2)
    f_end = zeros(d)
    f!(f_end, u, s, 1.)
    f_trj = [f_end f_trj[:,1:end-1]]

    y, xi = lss(du_trj, dJ_trj, 
    zeros(d,n), s, 2, f_trj)
    dJds = compute_sens(y, zeros(n), X_trj, zeros(d,n))
    vsh[:,:] = y
   
    save(string("../data/rijke_adjoint_sensitivity/",
        "vsh_and_dJds_", string(n_spe), ".jld"),
             "dJds", dJds, "vsh", vsh)
end
function collect_adjoint_sensitivities()
    pmap(Rijke_adjoint_sensitivity, 897:1120)
end


#=
    n_samples = 1000
    dJacds_ens, dJrayds_ens = zeros(2,n_samples), 
                    zeros(2,n_samples)
    up = zeros(d)
    um = zeros(d)
    eps = 1.e-3
    n = 100 # produces a growth of ≈ e
    for i = 1:n_samples
        @time dJacds_ens[:,i] = Zygote.gradient(s ->Jac(u, s, n), s)[1]
    dJrayds_ens[:,i] = Zygote.gradient(s ->Jray(u, s, n),
                               s)[1]
    u .= Rijke(u, [7.0, 0.2], n)
    end
=#
