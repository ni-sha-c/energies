include("lorenz63.jl")
include("clvs.jl")
"""
clvs : nxdxd_u
"""
	s = [10., 28., 8/3]
	m = 1
	n = 5000
	n_runup = 5000
	u0 = rand(3,m)
	u_init = lorenz63(u0, s, n)[end,:,:]
	u_trj = lorenz63(u_init, s, n)[:,:,1]
	du_trj = dlorenz63(u_trj, s)
	du_trj = permutedims(du_trj,[2,3,1])
	X = perturbation(u_trj,s) #ith col in T_{u_{i+1}} M
    
	d = size(du_trj)[1]
    n = size(du_trj)[3]
	d_u = 3

	lyap_exps = zeros(d_u)
    R = zeros(d_u,d_u,n)
    Q = zeros(d,d_u,n)
	v = zeros(d,1,n) #assume one parameter
	v[:,1,1] = rand(3)
	A = qr!(Q[:,:,1])
    Q[:,:,1] = Array(A.Q)
    R[:,:,1] = A.R
	
	b = zeros(d_u,1,n) #1 is the # parameters

    for i=2:n
		v[:,:,i] = du_trj[:,:,i-1]*v[:,:,i-1] + 
					X[:,i-1]
		Q[:,:,i] = du_trj[:,:,i-1]*Q[:,:,i-1]
        A = qr!(Q[:,:,i])
        Q[:,:,i] = Array(A.Q)
        R[:,:,i] = A.R
        b[:,:,i] = (Q[:,:,i]')*v[:,:,i]
		v[:,:,i] = v[:,:,i] - Q[:,:,i]*b[:,:,i]
		lyap_exps .+= log.(abs.(diag(R[:,:,i])))./n
    end

	# LSS solve
	using SparseArrays
	using BlockArrays
	b = reshape(collect(b),d_u,n)
	ndu = n*d_u
	eye = sparse(Matrix([zeros(ndu, d_u) 1.0I(ndu)]))
	D = (BlockArray{Float64}(zeros(ndu,ndu),
               d_u*ones(Int64,n),
               d_u*ones(Int64,n)))
	[D[Block(i,i)] = R[:,:,i] for i =1:n]
	D = sparse([Array(D) zeros(ndu, d_u)])
	B = D - eye
	BB = B*transpose(B)
	a = -transpose(B)*(BB\b[:])
	a = reshape(a, d_u, n+1)[:,1:end-1]
~
	# shadowing direction
	v = reshape(collect(v),d,n)
	vsh = zeros(d, n)
	for i = 1:n
		vsh[:,i] = v[:,i] + Q[:,:,i]*reshape(a[:,i],d_u,1)
	end

