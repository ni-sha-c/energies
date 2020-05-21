function lss(u_trj, du_trj, X, f, J, dJ, s, d_u)
	n, d = size(u_trj)
	lyap_exps = zeros(d_u)
    R = zeros(d_u,d_u,n)
    Q = zeros(d,d_u,n)
	v = zeros(d,1,n) #assume one parameter
	v[:,1,1] = rand(3)
	A = qr!(Q[:,:,1])
    Q[:,:,1] = Array(A.Q)
    R[:,:,1] = A.R
	
	b = zeros(d_u,1,n) #1 is the # parameters
	ff = zeros(n)
	[ff[i] = norm(f[:,i])^2.0 for i = 1:n]
    for i=2:n
		v[:,:,i] = du_trj[:,:,i-1]*v[:,:,i-1] + 
					X[:,i-1]
		Q[:,:,i] = du_trj[:,:,i-1]*Q[:,:,i-1]
		[Q[:,j,i] = Q[:,j,i-1] - (f[:,i]'*
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
	println("Solving the least squares problem...")
	include("lsssolve.jl")
	a = lsssolve(R,b)
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
	return vsh, dJds
end
