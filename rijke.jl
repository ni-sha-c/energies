# Rijke tube model
# Reference: Huhn and Magri 2020
using LinearAlgebra
include("cheb.jl")

function rijke(u0, beta, nSteps)
c1 = 0.1
c2 = 0.06
xf = 0.2
tau = 0.2
Ng = 10
tNg = 2*Ng + 1
Nc = 10
N = 2*Ng + Nc
j = LinRange(1,Ng,Ng)
jpi = @. pi.*j
cjpixf = @. cos.(xf.*jpi)
sjpixf = @. sin.(xf.*jpi)
zetaj = @. c1.*j.*j .+ c2.*
		sqrt.(j)
uffun(eta) = dot(eta,cjpixf)
qfun(t) = beta*(sqrt(abs(1.0 + t)) - 1.0)
#qfun(t) = beta*((1/3 + t)^2.0 + 1.e-4)^0.25 - 
#			1/sqrt(3)
D = cheb_diff_matrix(Nc-1)

function f(u::Array{Float64,1}, uDot::Array{Float64,1})
	eta = view(u, 1:Ng)
	mu = view(u, Ng+1:2*Ng)
	v = view(u, tNg:N)
	uDot[1:Ng] .= @. jpi.*mu
	uDot[Ng+1:2*Ng] .= @. -jpi.*eta - zetaj.*mu - 2.0.*
	qfun(v[1]).*sjpixf
	uDot[tNg:N] .= -2.0/tau.*(D*v)
end
# Time integration
# 3-stage scheme from Hager 2000 pg 272
#u0 = zeros(N)
#u0[1] = 1.
dt = 1.e-1/Nc*tau/2.0

uDot1 = zeros(N)
uDot2 = zeros(N)
uDot3 = zeros(N)

ufMean = zeros(nSteps)
@time for n = 1:nSteps
	ufMean[n] = norm(u0[1:2*Ng])^2.0/4.0
	#println("timestep ", n, " mean(u0) " , 
	#		uMean[n])
	u = copy(u0)
	f(u, uDot1)
	u .= u .+ (1/2)*dt.*uDot1
	f(u, uDot2)
	u .= u0 .- dt.*uDot1 .+ 
		2*dt.*uDot2
	f(u, uDot3)
	u0 .= u0 .+ dt.*((1/6).*uDot1 .+ 
					 (2/3).*uDot2 .+ (1/6).*uDot3)
	u0[end] = uffun(u0[1:Ng])
end
nRunUp = 20000
timeLength = nSteps - nRunUp + 1
return sum(ufMean[nRunUp:end])/timeLength
end
