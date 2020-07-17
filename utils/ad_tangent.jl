include("../examples/lorenz63.jl")
using Test
using DiffEqSensitivity, OrdinaryDiffEq, ForwardDiff
using Zygote
function time_average_obj(s, n, u0)
	u = copy(u0)
	J_avg = 0.
	for i = 1:n
		u = lorenz63_ad(u, s, 1)
		J_avg += u[end]/n
	end
	return J_avg
end
	s = [10., 28., 8/3]
	m = 1
    n = 200
    n_runup = 5000
	n_samples = 1000
    u0 = rand(3)
	u0 = lorenz63_ad(u0, s, n_runup)
	dJds = zeros(3,n_samples)
	for i = 1:n_samples
		u0 .= lorenz63_ad(u0, s, 1)
		dJds[:,i] .= ForwardDiff.gradient(s->time_average_obj(
										s,n,u0), s)
	end
	@show sum(dJds, dims=2)/n_samples

