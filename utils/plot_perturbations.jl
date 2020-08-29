include("../examples/rijke.jl")
using OrdinaryDiffEq
using DiffEqSensitivity, Zygote
using JLD
using Distributed
#using PyPlot
using SharedArrays
function ad_solve(u0, eps, n, v, ind)
	d = N
	s = [7.0, 0.2]
	un = u0 + eps*v
	un = Rijke_ODE(un, s, n)
	return un[ind]
end

function compute_perturbations()
d = N
u = rand(d)
s = [7.0, 0.2]
nSteps = 1000
ad_norm = zeros(nSteps)
ad = SharedArray{Float64}(d)
eps_ad = 0.
v_norm, w_norm, fd_norm = zeros(nSteps), 
				zeros(nSteps), zeros(nSteps)
v, w, fd = rand(d), rand(d),
				rand(d)
u_trj = zeros(d, nSteps)
eps = 1.e-4
u_p = u + eps*fd
for n = 1:nSteps
	ans = @distributed for j = 1:d
			ad[j] = Zygote.gradient(eps -> ad_solve(u, eps, 1, v, j), eps)[1]
	end
	wait(ans)
	ad_norm[n] = norm(ad)


	du = dRijke(u, s, 1.e-6)
	v = du*v
	v_norm[n] = norm(v)
	
	u_trj[:,n] = u
	u = Rijke_ODE(u, s, 1)
	u_p = Rijke_ODE(u_p, s, 1)

	fd = (u_p - u)./eps
	fd_norm[n] = norm(fd)


end
for n = nSteps:-1:1
	u = u_trj[:,n]
	du = dRijke(u, s, 1.e-6)
	du = du'
	w = du*w
	w_norm[n] = norm(w)
end
save("../data/rijke_perturbations/v_w_fd_ad_norms.jld", 
	 "v_norm", v_norm,
	 "w_norm", w_norm,
	 "fd_norm", fd_norm,
	 "ad_norm", ad_norm)
end
function plot_perturbations()
	X = load("../data/rijke_perturbations/v_w_fd_norms.jld")
	v_norm = X["v_norm"]
	w_norm = X["w_norm"]
	fd_norm = X["fd_norm"]
	X = load("../data/rijke_perturbations/ad_norms.jld")
	ad_norm = X["ad_norm"]
	nSteps = div(size(v_norm)[1],1)
	v_norm = v_norm[1:nSteps]
	fd_norm = fd_norm[1:nSteps]
	w_norm = w_norm[1:nSteps]
	ad_norm = ad_norm[1:nSteps]
	fig, ax = subplots(1,1)
	ax.semilogy(dt*(1:nSteps), v_norm, "v", ms=4.0,
				label="tangent")
	ax.semilogy(dt*(1:nSteps), w_norm, ".", ms=4.0,
				label="adjoint")
	ax.semilogy(dt*(1:nSteps), fd_norm, "^", ms=4.0,
				label="FD")
	#ax.semilogy(dt*(1:nSteps), ad_norm, "1", ms=5.0,
	#			label="AD")

	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax.set_xlabel("time", fontsize=28)
	ax.set_ylabel("Perturbation norms", fontsize=28)
	ax.grid(true)
	lgnd = fig.legend(loc=(0.3,0.6),
					  fontsize=28)
	lgnd.legendHandles[1]._legmarker.set_markersize(20)
	lgnd.legendHandles[2]._legmarker.set_markersize(20)
	lgnd.legendHandles[3]._legmarker.set_markersize(20)
	lgnd.legendHandles[4]._legmarker.set_markersize(20)
end

