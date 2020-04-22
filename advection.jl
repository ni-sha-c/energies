# Advection equation solution
# \partial v/\partial t + a \partial v/\partial x = 0
using PyPlot
using LinearAlgebra
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
for n = 1:nsteps
	plot(x,v1)
	pause(0.1)
	mul!(v1,A,v1)
end


