include("rijke.jl")
include("clvs.jl")
using PyPlot	
	s = [7.0, 0.2]
	#n = 100000
	n = 50000
	d = N
	nRunup = 40000
	u = Rijke(rand(d),s,nRunup)
	du_trj = zeros(d,d,n)
	duT_trj = zeros(d,d,n)
	f_trj = zeros(d,n)
	for i = 1:n
		f!(view(f_trj,:,i), u, s, 1.)
		du_trj[:,:,i] = dRijke(u, s, 1.e-8)
		duT_trj[:,:,n+1-i] = du_trj[:,:,i]'
		u .= Rijke(u, s, 1)
	end
	d_u = 20
	les, tan_clvs = clvs(du_trj, d_u, zeros(d,n))
	les1, adj_clvs = clvs(duT_trj, d_u, zeros(d,n))

	tan_angles = zeros(d_u, d_u, n)
	adj_angles = zeros(d_u, d_u, n)
	
	for i = 1:n
		tan_angles[:,:,i] = tan_clvs[:,:,i]'*tan_clvs[:,:,i]
		adj_angles[:,:,n+1-i] = adj_clvs[:,:,i]'*
								adj_clvs[:,:,i]
	end
	n_spinup = 5000
	tan_adj_angles = zeros(d_u, d_u, n-1)
	for i = 2:n
		tan_adj_angles[:,:,i-1] = tan_clvs[:,:,i]'*
		adj_clvs[:,:,n+2-i]
	end
	tan_adj_angles = tan_adj_angles[:,:,n_spinup:n-1-n_spinup]
	mean_tan_adj_angles = 180*acos.(sum(tan_adj_angles,
										dims=3)/n)/pi
	mean_tan_adj_angles = mean_tan_adj_angles[:,:,1]
								


	mean_tan_angles = sum(tan_angles,dims=3)/n
	mean_tan_angles = mean_tan_angles[1:6,1:6,1]
	mean_tan_angles = 180*acos.(mean_tan_angles)/pi
	mean_adj_angles = sum(adj_angles,dims=3)/n
	mean_adj_angles = mean_adj_angles[1:6,1:6,1]
	mean_adj_angles = 180*acos.(mean_adj_angles)/pi
	

	fig, ax = subplots(1,1)
	ax.xaxis.set_tick_params(labelsize=18)
	ax.yaxis.set_tick_params(labelsize=18)
	angle_plot = ax.matshow(mean_tan_angles)
	cbar = fig.colorbar(angle_plot)
	cbar.ax.tick_params(labelsize=18)
	fig2, ax2 = subplots(1,1)
	ax2.xaxis.set_tick_params(labelsize=18)
	ax2.yaxis.set_tick_params(labelsize=18)
	angle_plot2 = ax2.matshow(mean_adj_angles)
	cbar = fig2.colorbar(angle_plot2)
	cbar.ax.tick_params(labelsize=18)


	fig1, ax1 = subplots(1,1)
	ax1.xaxis.set_tick_params(labelsize=18)
	ax1.yaxis.set_tick_params(labelsize=18)
	ax1.set_title(L"Probability distribution of 
				  angles between $V^1$ and $V^2$", 
				  fontsize=18)
	ax1.hist(acos.(tan_angles[1,2,:])*180/pi,
			 density=true,bins=100)
	# minimum angle between 1 and 2 is 15.6 degrees
	ax1.set_xlabel("Angle in degrees",fontsize=18)
	ax1.grid(true)
	fig3, ax3 = subplots(1,1)
	ax3.xaxis.set_tick_params(labelsize=18)
	ax3.yaxis.set_tick_params(labelsize=18)
	ax3.set_title(L"Probability distribution of 
				  angles between $W^1$ and $W^2$", 
				  fontsize=18)
	ax3.hist(acos.(adj_angles[1,2,:])*180/pi,
			 density=true,bins=100)
	# minimum angle between 1 and 2 is 4 degrees
	ax3.set_xlabel("Angle in degrees",fontsize=18)
	ax3.grid(true)

	fig4, ax4 = subplots(1,1)
	ax4.xaxis.set_tick_params(labelsize=18)
	ax4.yaxis.set_tick_params(labelsize=18)
	angle_plot4 = ax4.matshow(mean_tan_adj_angles)
	cbar = fig4.colorbar(angle_plot4)
	cbar.ax.tick_params(labelsize=18)

