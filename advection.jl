# Advection equation solution
# \partial v/\partial t + a \partial v/\partial x = 0
using PyPlot
using LinearAlgebra
include("cheb.jl")
nsteps = 100
tau = 0.02
a = 2.0/tau
dx = 0.05
dt = dx/a*0.9
m = Int(cld(2,dx))
x = LinRange(-1.,1.,m)
v0 = sin.(pi*x)
A = Bidiagonal(ones(m), -ones(m-1), :L)
c = dt*a/dx
A .= lmul!(c,A)
A .= I(m) .- A  
v1 = copy(v0)
Nc = 20
D = cheb_diff_matrix(Nc)
dtc = 0.1/Nc/a
D .= I(Nc+1) - a.*dtc.*D
xc = cheb_pts(Nc)
p0 = sin.(pi*xc)
p1 = copy(p0)
fig, ax = subplots(1,2)
for n = 1:1000
	ax[1].plot(x,v1)
	ax[2].plot(xc,p1)
	pause(0.1)
	mul!(v1,A,v1)
	p1 .= D*p1
end


