include("rijke.jl")
u = zeros(N)
u[1] = 1.
s = [2.5, 0.2]
u_trj, uf_trj = Rijke(u, s, 200000) 
#print(u_trj[1:10])
#sol = Rijke_ODE(u, s, 1)
