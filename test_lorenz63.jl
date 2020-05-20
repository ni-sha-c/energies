include("lorenz63.jl")
function test_dlorenz63()
	s = [10., 28., 8/3]
	m = 10
	u0 = rand(3,m)
	u_trj = lorenz63(u0, s, 1)[2,:,:]
	eps = 1.e-4
	du_x = (lorenz63(u_trj .+ eps.*[1; 0.; 0.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[1; 0.; 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_y = (lorenz63(u_trj .+ eps.*[0; 1.; 0.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[0; 1.; 0.], 
						s, 1)[2,:,:])./(2*eps)
	du_z = (lorenz63(u_trj .+ eps.*[0; 0.; 1.],
						s, 1)[2,:,:] - 
				lorenz63(u_trj .- eps.*[0; 0.; 1.], 
						s, 1)[2,:,:])./(2*eps)
	du_fd = reshape(collect([du_x' du_y' du_z']), 
						m, 3, 3)
	du = dlorenz63(u_trj', s)
	@assert all(isapprox.(du_fd, du)) == 1
end

