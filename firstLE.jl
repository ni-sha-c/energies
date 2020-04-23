include("rijke.jl")
N = 30
u0Init = zeros(N)
u0Init[1] = 1.
beta = 7.0
nRunUp = 500000
u0 = rijke(u0Init, beta, nRunUp)
n_samples = 50000
u1 = zeros(N)
v0 = rand(N)
v0 ./= norm(v0)
v1 = zeros(N)
eps = 1.e-2
nSeg = 5
l1 = zeros(n_samples)
#dt = 1.e-1/10/2.0
dt = 1
for i = 1:n_samples
	u0Orig = copy(u0)
	u1 .= rijke(u0 .+ eps*v0, beta, nSeg)
	v1 .= (u1 .- u0Orig)./eps
	u0 .= u1
	l1[i] = log(norm(v1))/nSeg/dt
end
#plot(cumsum(l1)./StepRange(1,1,n_samples))

