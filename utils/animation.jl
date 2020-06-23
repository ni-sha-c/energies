# initialize the attractor
n = 100
n = 10000
s = [10., 28., 8/3]
u0 = rand(3,1)
u_trj = lorenz63(u0, s, n)
# initialize a 3D plot with 1 empty series
plt = path3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
                xlab = "x", ylab = "y", zlab = "z",
                title = "Lorenz Attractor", marker = 1)

# build an animated gif, saving every 10th frame
u_trj = permutedims(u_trj,[3,2,1])
@gif for i=1:n
		push!(plt, u_trj[i,1,:], u_trj[i,2,:], u_trj[i,3,:])
end every 10
