# LSS solve
using SparseArrays
using BlockArrays
function lsssolve(R, b, vf, Qf)
    d_u, n = size(b)
    ndu = n*d_u
    eye = sparse(Matrix([zeros(ndu, d_u) 1.0I(ndu)]))
    D = (BlockArray{Float64}(zeros(ndu,ndu),
               d_u*ones(Int64,n),
               d_u*ones(Int64,n)))
    [D[Block(i,i)] = R[:,:,i] for i =1:n]
    D = sparse([Array(D) zeros(ndu, d_u)])
    B = D - eye
	pf_Q = [Qf zeros(d_u,1)][:]
	B = sparse([Array(B); pf_Q'])
    BB = B*transpose(B)
	b1 = [b[:]; vf]
	a = -transpose(B)*(BB\b1[:])
    a = reshape(a, d_u, n+1)[:,1:end-1]
	return a
end
