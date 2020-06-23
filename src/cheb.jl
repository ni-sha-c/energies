# chebyshev points and chebyshev differentiation matrix 
# Ref: Spectral methods in MATLAB by Trefethen, chapters 1, 5 and 6.
function cheb_pts(n)
	return cos.((0:n)*pi/n)
end
function cheb_diff_matrix(n)
	x = cheb_pts(n)
	np = n+1
	D = zeros(np,np)
	E = view(D, reverse(1:np*np))
	D[1] = (2*n*n + 1)/6.
	E[1] = -D[1]


	nhalf = fld(np, 2)
	xint = view(x,2:n)
	D[2:n] .= (-1).^(1:n-1)./(xint .- 1)./2.0
	E[2:n] .= -D[2:n]
	D[np] = (-1)^np/2.0
	E[np] = -D[np]

	if mod(np,2) != 0
		j = nhalf + 1
		A = view(D, (j-1)*np+1:j*np)
		A .=  1.0./(x .- x[j]).*
		((-1).^((1 + j):(np + j)))
		A[j] = -x[j]/2/(1. - x[j]*x[j])
		A[1] *= 2.0
		A[np] *= 2.0
	end
	
	for j=2:nhalf
		A = view(D, (j-1)*np+1:j*np)
		B = view(E, (j-1)*np+1:j*np)
		A .= 1.0./(x .- x[j]).*
		((-1).^((1 + j):(np + j)))
		B .= -A
		A[1] *= 2.0
		B[1] *= 2.0
		A[j] = -x[j]/2/(1. - x[j]*x[j])
		B[j] = -A[j]
		A[np] *= 2.0
		B[np] *= 2.0
	end
	
	return D
end
