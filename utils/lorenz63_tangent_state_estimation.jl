include("../examples/lorenz63.jl")
include("../src/lss.jl")
#create epsilon orbit.
u = rand(3,1)
n_runup = 5000
#u = lorenz63(u, 
