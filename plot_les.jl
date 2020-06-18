include("rijke.jl")
include("clvs.jl")
d = 30
u = rand(d)
nRunup = 40000
s = [7.0,0.2]
u = Rijke_ODE(u, s, nRunup)[:,end]
println("Done with Runup")
nSteps = 1
u0 = copy(u)
du = zeros(d,d)
function one_step(u0, s=[7.0,0.2])
	return Rijke_ODE(u0, s, 1)
end
@show u0
for i = 1:nSteps
@time for j = 1:d
		a = Zygote.gradient(v0->one_step(u0 + v0)[j,end],u0) 
		du[j,:] = a[1]
	end
end
@show u0
#println("Done with computing primal trajectory")
# get Jacobian.
@time dTu = dRijke(u0, s, 1.e-8)
#println("Done with getting Jacobian trajectory")
#println("Getting CLVs...")
# get clvs
#les, clv_trj = clvs(dTu, 20)

