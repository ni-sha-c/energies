include("../examples/rijke.jl")
include("../src/clvs.jl")
d = 30
u = rand(d)
nRunup = 1000000
s = [7.0,0.2]
u = Rijke_ODE(u, s, nRunup)
println("Done with Runup")
nSteps = 2000
du_trj = zeros(d,d,nSteps)
du_trj_ad = zeros(d,d,nSteps)
delayed_velocity = zeros(nSteps)
dheat = zeros(nSteps)
for i = 1:nSteps
	u .= Rijke_ODE(u, s, 1)
	delayed_velocity[i] = u[2*Ng + 1]
	du_trj[:,:,i] = dRijke(u,s,1.e-5)
end
function dqfun(t)
	if abs(t + 1.0) > 0.01
		dqdot_dt = 0.5/sqrt(abs(1. + t))
		if 1. + t < 0.
			dqdot_dt *= -1.
		end
		return dqdot_dt
	end
	dqdot_dt = 0.
	for i = 1:5
			dqdot_dt += coeffs[i]*(i-1)*((1+t)^(i-2))
	end
	return dqdot_dt
end
#=
function du_test1(u, s, eps=0.)
	du_ad = zeros(d,d)
	for j = 1:d
		v = zeros(d)
		du_ad[:,j] = Zygote.gradient(v->
				Rijke_ODE(u .+ v, s, 1)[j],v)[1]
	end
	return du_ad'
end
=#
println("Done with computing primal and Jacobian 
		trajectory")
println("Getting CLVs...")
# get clvs
les, clv_trj = clvs(du_trj, 10)

