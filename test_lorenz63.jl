include("lorenz63.jl")
s = [10., 28., 8/3]
n = 1
m = 10
u0 = rand(3,m)
u_trj = lorenz63(u0, s, n)[:,:,1]
n = n + 1
u_p_x = u_trj .+ eps*ones(1,3)
u_m_x = u_trj
