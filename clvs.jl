using LinearAlgebra
"""
    Compute Lyapunov Exponents and Covariant 
    Lyapunov vectors.
    Reference: Ginelli 2013
	Inputs:
		DTu: dxdxm m-length timeseries of jacobian matrices
		du: tangent subspace dimension that the computed CLVs must span.
	Outputs:
		lyap_exps: the first du Lyapunov exponents
		Q: dxduxm m-length timeseries of the first du CLVs

"""
function clvs(DTu::Array{Float64,3},du::Int64,F::Array{Float64,2})
    d = size(DTu)[1]
    m = size(DTu)[3]
    lyap_exps = zeros(du)
    R = zeros(du,du,m)
    Q = zeros(d,du,m)
	println(size(Q), size(R))
    A = qr!(rand(d,du))
	Q[:,:,1] = Array(A.Q)
	R[:,:,1] = A.R
	FF = ones(m)
	if isapprox.(sum(F,dims=1),0.) == false
		FF = [sum(x->x*x, F[:,i]) for i=1:m]
	end
    for i=2:m
        Q[:,:,i] = DTu[:,:,i-1]*Q[:,:,i-1]
		[Q[:,j,i] = Q[:,j,i] .- dot(Q[:,j,i], 
					F[:,i])*F[:,i]/FF[i] for j=1:du]
        A = qr!(Q[:,:,i])
		Q[:,:,i] = Array(A.Q)
		R[:,:,i] = A.R
        lyap_exps .+= log.(abs.(diag(R[:,:,i])))./m
    end
    C = zeros(du,du,m)
    C[:,:,end] = diagm(ones(du))
    for i=reverse(1:m-1)
        C[:,:,i] = R[:,:,i+1]\C[:,:,i+1]
        [normalize!(view(C,:,j,i)) for j = 1:du]
        Q[:,:,i] = Q[:,:,i]*C[:,:,i]
    end
    return lyap_exps, Q
end

