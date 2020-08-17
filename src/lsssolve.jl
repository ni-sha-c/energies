# LSS solve
using SparseArrays
using BlockArrays
function lsssolve(R, b)
    d_u, n = size(b)
    ndu = n*d_u
    eye = sparse(Matrix([zeros(ndu, d_u) 1.0I(ndu)]))
    D = (BlockArray{Float64}(zeros(ndu,ndu),
               d_u*ones(Int64,n),
               d_u*ones(Int64,n)))
    [D[Block(i,i)] = R[:,:,i+1] for i =1:n-1]
    D = sparse([Array(D) zeros(ndu, d_u)])
    B = D - eye
    BB = B*transpose(B)
	a = -transpose(B)*(BB\b[:])
    a = reshape(a, d_u, n+1)[:,1:end-1]
	return a
end
