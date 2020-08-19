include("adjoint_lsssolve.jl")
"""
Implements the NILSAS algorithm due to Ni and Talnikar 2019.
"Adjoint sensitivity analysis on chaotic dynamical
systems by Non-Intrusive Least Squares Adjoint
Shadowing (NILSAS)"

Below, u -> F(u) is the time-1 primal map and u_trj refers 
to a primal trajectory of length n.
Inputs:
    du_trj: dxdxn n-length timeseries of Jacobian matrices
    		along u_trj
    X: dxn n-length timeseries of perturbation vectors. 
    		X[:,i] = dF/ds(u_trj[i,:]) which is in the 
    		tangent space of u_trj[i+1,:].
    f: dxn n-length timeseries of the neutral CLV. By 
    		assumption, the only neutral CLV is the RHS of 
    		the primal ode: du/dt = f(u,s). 
    		In case of maps, u_n = F(u_{n+1}), f is zeros.
    J:mxn n-length timeseries of m scalar functions J(u)
    		evaluated along u_trj.
    dJ: mxdxn n-length timeseries of DJ(u) evaluated along 
    		u_trj
    s: array of parameters
    d_u: integer dimension of the unstable subspace. You can 		set a guess value to obtain correct sensitivities 
    	as long as the guess value is greater than or equal
    	to the true d_u.
Outputs:
    vsh: dxn n-length timeseries of shadowing tangent 
    	vectors that are bounded solutions of the 
    	conventional tangent equation: 
    		v[:,i+1] = du_trj[:,i]*v[:,i] + X[:,i]
    
    dJds: scalar-valued sensitivity of interest:
    	d/ds sum(J) ≈ d/ds <J, mu>, where mu is 
    	the SRB measure.
"""
function lss(du_trj, X, f, s, d_u, g)
    d, n = size(X)
    lyap_exps = zeros(d_u)
    R = zeros(d_u,d_u,n)
    Q = zeros(d,d_u,n)
    v = zeros(d,1,n) #assume one parameter
    v[:,1,1] = zeros(d)
    A = qr!(rand(d,d_u))
    Q[:,:,1] = Array(A.Q)
    R[:,:,1] = A.R
    
    b = zeros(d_u,1,n) #1 is the # parameters
    pf_Q = zeros(d_u,1,n) #1 is the # parameters
    pf_v = zeros(n)

    ff = zeros(n)
    if all(isapprox.(f,0.)) 
    	ff = ones(n)
    else
    	[ff[i] = sum(x -> x^2, f[:,i]) for i = 1:n]
    end
    pf_Q[:,1,1] = Q[:,:,1]'*g[:,1]
    [Q[:,j,1] = Q[:,j,1] - (pf_Q[j,1,1]*
     		f[:,1]/ff[1]) for j=1:d_u]

    for i=2:n
        v[:,:,i] = du_trj[:,:,i-1]*v[:,:,i-1] + 
    				X[:,i-1]
    	Q[:,:,i] = du_trj[:,:,i-1]*Q[:,:,i-1]
    	pf_Q[:,:,i] = Q[:,:,i]'*g[:,i]
    	#[Q[:,j,i] = Q[:,j,i] - pf_Q[j,1,i]*f[:,i]/
    	#				ff[i] for j=1:d_u]
    	A = qr!(Q[:,:,i])
        Q[:,:,i] = Array(A.Q)
        R[:,:,i] = A.R
    	pf_v[i] = dot(v[:,1,i],g[:,i])
    	v[:,:,i] = v[:,:,i] - (pf_v[i]/ff[i]*f[:,i])
        b[:,:,i] = (Q[:,:,i]')*v[:,:,i]
    	v[:,:,i] = v[:,:,i] - Q[:,:,i]*b[:,:,i]
	lyap_exps .+= log.(abs.(diag(R[:,:,i])))./n
    end
    
    b = reshape(collect(b),d_u, n)
    a = lsssolve(R,b,sum(pf_v),pf_Q[:,1,:])
    v = reshape(collect(v),d,n)
    vsh = zeros(d, n)
    xi = zeros(n)
    vsh[:,1] = Q[:,:,1]*reshape(a[:,1], d_u, 1) + 
    			  v[:,1]
    for i = 2:n
    	vsh[:,i] = v[:,i] + Q[:,:,i]*reshape(a[:,i],d_u,1)  
    end
    return vsh, xi
end
function compute_sens(vsh, xi, dJ, f)
    m, d, n = size(dJ)
    dJds = zeros(m)
    for i = 1:n
    	dJds .= dJds .+ dJ[:,:,i]*vsh[:,i] 
    end
    return dJds
end
