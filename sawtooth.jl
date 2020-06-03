function sawtooth(u0, s=0., n=1)
	m, = size(u0)
	n = n+1
	u_trj = zeros((m,n))
	u_trj[:,1] = u0
	for i = 2:n
		x = view(u_trj,:,i-1)
		u_trj[:,i] = (2*x + s*sin.(2*pi*x/4)) .% 1
	end 
	return permutedims(u_trj,[2,1])
end
function dsawtooth(u, s=0)
	n, = size(u)
	return 2.0 .+ s*cos.(2*pi*u/4)*2*pi/4
end
function perturbation(u,s)
	n, d = size(u)
	# the perturbation in row i in T_{u_(i+1)} M
	return [zeros(1,n); dt*u[:,1]'; zeros(1,n)]
end

