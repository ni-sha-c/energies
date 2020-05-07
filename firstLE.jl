include("rijke.jl")
#using Plots
N = 30
u0Init = zeros(N)
u0Init[1] = 1.
s = [7.0, 0.2]
nRunUp = 40000
u0 = Rijke(u0Init, s, nRunUp)
n_samples = 8
u1 = zeros(N)
v0 = rand(N)
v0 ./= norm(v0)
v1 = zeros(N)
eps = 1.e-4
l1 = zeros(n_samples)
#dt = 1.e-1/10/2.0
dt = 1.e-3
nSeg = 1
for i = 1:n_samples
	u0Init .= copy(u0)
	u1 .= Rijke(u0 .+ eps*v0, s, nSeg)
	v1 .= (u1 .- u0Init)./eps
	u0 .= u1
	l1[i] = log(norm(v1))/nSeg/dt
end
#plot(cumsum(l1)./StepRange(1,1,n_samples))

