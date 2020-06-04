include("lsssolve.jl")
"""
Implements the NILSS algorithm due to Ni and Wang JCP 2017
Below, u -> F(u) is the time-1 primal map.
Inputs:
	u_trj: nxd n-length timeseries of primal solutions
	du_trj: dxdxn n-length timeseries of Jacobian matrices
			along u_trj
	X: dxn n-length timeseries of perturbation vectors. 
			X[:,i] = dF/ds(u_trj[i,:]) which is in the 
			tangent space of u_trj[i+1,:].
	f: dxn n-length timeseries of the neutral CLV. By 
			assumption, the only neutral CLV is the RHS of 
			the primal ode: du/dt = f(u,s). 
			In case of maps, u_n = F(u_{n+1}), f is zeros.
	J: n-length timeseries of the scalar function J(u)
			evaluated along u_trj
	dJ: dxn n-length timeseries of DJ(u) evaluated along 
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
function lss(u_trj, du_trj, X, f, J, dJ, s, d_u)
	n, d = size(u_trj)
	lyap_exps = zeros(d_u)
    R = zeros(d_u,d_u,n)
    Q = zeros(d,d_u,n)
	v = zeros(d,1,n) #assume one parameter
	v[:,1,1] = zeros(d)
	A = qr!(rand(d,d_u))
    Q[:,:,1] = Array(A.Q)
    R[:,:,1] = A.R
	
	b = zeros(d_u,1,n) #1 is the # parameters
	ff = zeros(n)
	if all(isapprox.(f,0.)) 
		ff = ones(n)
	else
		[ff[i] = sum(x -> x^2, f[:,i]) for i = 1:n]
	end
   	
	for i=2:n
	    v[:,:,i] = du_trj[:,:,i-1]*v[:,:,i-1] + 
					X[:,i-1]
		Q[:,:,i] = du_trj[:,:,i-1]*Q[:,:,i-1]
		[Q[:,j,i] = Q[:,j,i] - (f[:,i]'*
						Q[:,j,i])*f[:,i]/
		 				ff[i] for j=1:d_u]
		A = qr!(Q[:,:,i])
        Q[:,:,i] = Array(A.Q)
        R[:,:,i] = A.R
        b[:,:,i] = (Q[:,:,i]')*v[:,:,i]
		v[:,:,i] = (v[:,:,i] - 
					Q[:,:,i]*b[:,:,i] - 
					(f[:,i]'*v[:,:,i])/
					(ff[i]).*f[:,i])
		lyap_exps .+= log.(abs.(diag(R[:,:,i])))./n
    end
	
	
	b = reshape(collect(b),d_u, n)
	println("Solving the least squares problem... ")
	a, condno = lsssolve(R,b)
~
	# shadowing direction
	println("Computing shadowing direction...")
	v = reshape(collect(v),d,n)
	vsh = zeros(d, n)
	xi = zeros(n)
	for i = 1:n
		vsh[:,i] = v[:,i] + Q[:,:,i]*reshape(a[:,i],d_u,1)
		xi[i] = vsh[:,i]'*f[:,i]/ff[i]
	end

	# sensitivity
	println("Computing sensitivity...")
	dJds = 0.
	Jmean = sum(J)/n
	for i = 1:n
			dJds += vsh[:,i]'*dJ[:,i]/n + 
			xi[i]*(Jmean .- J[i])
	end
	return vsh, dJds, condno
end

