include("rijke.jl")
include("clvs.jl")
d = 30
u = rand(d)
nRunup = 40000
s = [7.0,0.2]
u = Rijke(u, s, nRunup)
println("Done with Runup")
nSteps = 100000
u_trj = zeros(d,nSteps)
u_trj[:,1] = u
for i = 2:nSteps
	u_trj[:,i] = Rijke(u_trj[:,i-1], 
				s, 1)
end
println("Done with computing primal trajectory")
# get Jacobian.
@time dTu = dRijke(u_trj, s, 1.e-4)
println("Done with getting Jacobian trajectory")
println("Getting CLVs...")
# get clvs
les, clv_trj = clvs(dTu, 20)

