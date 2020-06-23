using PyPlot
include("rijke.jl")
fig, ax = subplots(2,2)
dt = 1/40.
beg = 800
endd = 860
tstart = Int(beg/dt)
tend = Int(endd/dt)
timeStepArr = StepRange(tstart, 10, tend) 
t = LinRange(0, endd-beg, length(timeStepArr)) 
println("Doing beta = ", 0.4)
uf1 = rijke(0.4)

println("Doing beta = ", 2.5)
uf2 = rijke(2.5)

println("Doing beta = ", 7.0)
uf3 = rijke(7.0)

println("Doing beta = ", 10.0)
uf4 = rijke(10.0)

ax[1].set_title(L"$\beta = 0.4$",fontsize=20)
ax[2].set_title(L"$\beta = 2.5$",fontsize=20)
ax[3].set_title(L"$\beta = 7.0$",fontsize=20)
ax[4].set_title(L"$\beta = 10.0$",fontsize=20)

ax[1].plot(t, uf1[tstart:10:tend])
ax[2].plot(t, uf2[tstart:10:tend])
ax[3].plot(t, uf3[tstart:10:tend])
ax[4].plot(t, uf4[tstart:10:tend])

ax[1].xaxis.set_tick_params(labelsize=20)
ax[2].xaxis.set_tick_params(labelsize=20)
ax[3].xaxis.set_tick_params(labelsize=20)
ax[4].xaxis.set_tick_params(labelsize=20)

ax[1].yaxis.set_tick_params(labelsize=20)
ax[2].yaxis.set_tick_params(labelsize=20)
ax[3].yaxis.set_tick_params(labelsize=20)
ax[4].yaxis.set_tick_params(labelsize=20)

ax[2].set_xlabel(L"$t$",fontsize=20)
ax[4].set_xlabel(L"$t$",fontsize=20)
ax[1].set_ylabel(L"$u_f$",fontsize=20)
ax[2].set_ylabel(L"$u_f$",fontsize=20)
