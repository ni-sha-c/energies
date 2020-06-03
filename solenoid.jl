dt = 0.005
function solenoid(u0, s, n)
	s0, s1 = s
	m, d = size(u0)
	n = n+1
	u_trj = zeros((m,d,n))
	u_trj[:,:,1] = u0
	for i = 2:n
		x = view(u_trj,:,1,i-1)
		y = view(u_trj,:,2,i-1)
		z = view(u_trj,:,3,i-1)
	
		r, t = cart_to_cyl(x,y)

		r_next = @. s0 + (r - s0)/s1 + cos(t)/2
		t_next = @. 2*t
		z_next = @. z/s1 + sin(t)/2

		u_trj[:,1,i], u_trj[:,2,i] = cyl_to_cart(r_next,
												 t_next)
		u_trj[:,3,i] = z_next
	end 
	return permutedims(u_trj,[3,2,1])
end
function dsolenoid(u, s)
	s0, s1 = s
	n, d = size(u)
	x = view(u,:,1)
	y = view(u,:,2)
	z = view(u,:,3)
	
	T1r, T1t = cart_to_cyl(x,y)
	dT1 = dcart_to_cyl(x,y)
	dT1_dx = view(dT1, :, :, 1)
	dT1_dy = view(dT1, :, :, 2)
	
	T2r, T2t = s0 .+ (T1r .- s0)/s1 .+ cos.(T1t)/2, 
				2*T1t
	dT2r_dT1 = [1/s1*ones(n) -sin.(T1t)/2]
	dT2t_dT1t = 2.0*ones(n)
	dT2z_dT1t = cos.(T1t)/2 
	dT2z_dz = 1/s1*ones(n)
	
	dT3 = dcyl_to_cart(T2r, T2t)
	dT3_dT2r, dT3_dT2t = dT3[:,:,1], dT3[:,:,2]
	
	dT2r_dx = sum(dT2r_dT1.*dT1_dx, dims=2)
	dT2t_dx = dT2t_dT1t.*dT1_dx[:,2]
	dT2z_dx = dT2z_dT1t.*dT1_dx[:,2]

	dT2r_dy = sum(dT2r_dT1.*dT1_dy, dims=2)
	dT2t_dy = dT2t_dT1t.*dT1_dy[:,2]
	dT2z_dy = dT2z_dT1t.*dT1_dy[:,2]

	dT3_dx = [dT3_dT2r zeros(n)].*
			 [dT2r_dx dT2r_dx zeros(n)] .+ 
			 [dT3_dT2t zeros(n)].*
			 [dT2t_dx dT2t_dx zeros(n)] .+ 
			 [zeros(n) zeros(n) dT2z_dx]
	dT3_dy = [dT3_dT2r zeros(n)].*
			 [dT2r_dy dT2r_dy zeros(n)] .+ 
			 [dT3_dT2t zeros(n)].*
			 [dT2t_dy dT2t_dy zeros(n)] .+ 
			 [zeros(n) zeros(n) dT2z_dy]
	dT3_dz = [zeros(n) zeros(n) dT2z_dz]
			 
	return reshape([dT3_dx'; dT3_dy'; dT3_dz'],
				   d,d,:)
end
function perturbation(u,s)
	n, d = size(u)
	s0, s1 = s
	# the perturbation in row i in T_{u_(i+1)} M
	x, y, z = view(u,:,1), view(u,:,2), view(u,:,3)
	r, t = cart_to_cyl(x,y)
	r1, t1 = s0 .+ (r .- s0)/s1 + cos.(t)/2,
			 2*t
	# dTu_dr1 must be dim 2xn		 
	dTu_dr1 = transpose(dcyl_to_cart(r1, t1)[:,:,1]) 
	dr1_ds = 1.0 - 1/s1
	return [dr1_ds*dTu_dr1; zeros(1,n)]
end
function cart_to_cyl(x, y)
	return [sqrt.(x.*x .+ y.*y), 
			mod.(pi/180*atand.(y,x) .+ 2*pi, 2*pi)]
end
function cyl_to_cart(r,t)
	return [r.*cos.(t), r.*sin.(t)]
end
function dcart_to_cyl(x, y)
	r2 = x.*x + y.*y 
	r = sqrt.(r2)
	return reshape([x./r -y./r2 y./r x./r2],:,2,2)
end
function dcyl_to_cart(r, t)
	return reshape([cos.(t) sin.(t) -r.*sin.(t) r.*cos.(t)]
				   ,:, 2, 2)
end
function objective(u,s)
	n, d = size(u)
	u = u'
	J = zeros(n)
	[J[i] = sqrt(u[1,i]*u[1,i] + 
				 u[2,i]*u[2,i] + 
				 u[3,i]*u[3,i]) for i=1:n]
	return J
end
function dobjective(u,s)
	n, d = size(u)
	x, y, z = u[:,1], u[:,2], u[:,3]
	r = sqrt.(x.*x + y.*y + z.*z)
	dJ = reshape([x./r y./r z./r],n,d)
	dJ = dJ'
	return dJ
end
