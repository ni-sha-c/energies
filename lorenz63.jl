function lorenz63(u0, s, n)
	sigma, rho, beta = s
	dt = 0.005
	d, m = size(u0)
	n = n+1
	u_trj = zeros((m,d,n))
	u_trj[:,:,1] = u0'
	for i = 2:n
		x = u_trj[:,1,i-1]
		y = u_trj[:,2,i-1]
		z = u_trj[:,3,i-1]

		u_trj[:,1,i] = x + dt*(sigma.*(y - x))
		u_trj[:,2,i] = y + dt*(x.*(rho .- z) - y)
		u_trj[:,3,i] = z + dt*(x.*y - beta.*z)
	end 
	return permutedims(u_trj,[3,2,1])
end
function dlorenz63(u, s)
	sigma, rho, beta = s
	n, d = size(u)
	x = view(u,:,1)
	y = view(u,:,2)
	z = view(u,:,3)
	du = zeros(n, d, d)
	dt = 0.005
	@. du[:,1,1] = 1.0 - dt*sigma
	@. du[:,1,2] = dt*sigma
	@. du[:,2,1] = dt*(rho - z) 
	@. du[:,2,2] = 1.0 - dt
	@. du[:,2,3] = -dt*x 
	@. du[:,3,1] = dt*y
	@. du[:,3,2] = dt*x
	@. du[:,3,3] = 1.0 - dt*beta
	return du
end
#=
    def objective(self, fields, parameter):
        return fields[-1]

    def source(self, fields, parameter):
        sourceTerms = np.zeros_like(fields)
        sourceTerms[1] = self.dt*fields[0]
        return sourceTerms
        
    def gradientObjective(self, fields, parameter):
        dJ = np.zeros_like(fields)
        dJ[-1] = 1.0
        return dJ

    def tangentSolver(self, initFields, initPrimalFields, \
            parameter, nSteps, homogeneous=False):
        primalTrj = np.empty(shape=(nSteps, initFields.shape[0]))
        objectiveTrj = np.empty(nSteps)
        dt = self.dt
        primalTrj[0] = initPrimalFields
        objectiveTrj[0] = self.objective(primalTrj[0],parameter)
        for i in range(1, nSteps):
            primalTrj[i], objectiveTrj[i] = self.primalSolver(\
                    primalTrj[i-1], parameter, 1)
        xt, yt, zt = initFields
        sensitivity = np.dot(initFields, \
                    self.gradientObjective(primalTrj[0], parameter))/nSteps

        for i in range(nSteps):
            x, y, z = primalTrj[i]
            dxt_dt = self.sigma*(yt - xt) 
            dyt_dt = (parameter + self.rho - z)*xt - zt*x - yt 
            dzt_dt = x*yt + y*xt - self.beta*zt
            
            xt += dt*dxt_dt 
            yt += dt*dyt_dt
            zt += dt*dzt_dt
            
            finalFields = np.array([xt, yt, zt])
            if(homogeneous==False):
                finalFields += self.source(primalTrj[i],\
                        parameter)
                xt, yt, zt = finalFields

            if(i < nSteps-1):
                sensitivity += np.dot(finalFields, \
                    self.gradientObjective(primalTrj[i+1], parameter))/nSteps
        return finalFields, sensitivity
            

    def adjointSolver(self, initFields, initPrimalFields, \
            parameter, nSteps, homogeneous=False):
        rho = self.rho
        beta = self.beta
        sigma = self.sigma
        dt = self.dt
        primalTrj = np.empty(shape=(nSteps, initFields.shape[0]))
        objectiveTrj = np.empty(nSteps)

        primalTrj[0] = initPrimalFields
        objectiveTrj[0] = self.objective(primalTrj[0],parameter)
        for i in range(1, nSteps):
            primalTrj[i], objectiveTrj[i] = self.primalSolver(\
                    primalTrj[i-1], parameter, 1)
        xa, ya, za = initFields
        sensitivity = 0.
        for i in range(nSteps-1, -1, -1):
            x, y, z = primalTrj[i]
            dxa_dt = -sigma*xa + (parameter + rho - z)*ya + \
                    y*za 
            dya_dt = sigma*xa - ya + x*za 
            dza_dt = -x*ya - beta*za 
            
            xa += dt*dxa_dt 
            ya += dt*dya_dt
            za += dt*dza_dt
           
            finalFields = np.array([xa, ya, za])
            if(homogeneous==False):
                finalFields += self.gradientObjective(primalTrj[i],\
                        parameter)/nSteps
                xa, ya, za = finalFields
            if(i > 0):
                sensitivity += np.dot(finalFields, self.source(\
                    primalTrj[i-1], parameter))
        return finalFields, sensitivity
=#            



            
