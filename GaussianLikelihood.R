rm(list=ls(all=TRUE))
graphics.off()

GaussianLikelihood <- function(data,u,sigma) {

    n = length(data)
    par(mfrow = c(10,10),mar =  c(.5,.5,.5,.5))
    grid= 10
    
    sigmaf = sigma[1]
    sigmat = sigma[2]
    sigma = rev(seq(from = sigmaf, to = sigmat, by = (sigmat-sigmaf)/(grid-1)))

    meanf = u[1]
    meant = u[2]
    u = seq(from = meanf, to = meant, by = (meant-meanf)/(grid-1))

    gaussf = 0
    gausst = 2
    gauss = seq(from = gaussf, to = gausst, by = (gausst- gaussf)/90)
    Px_g_u_s = array(list(), dim = c(grid,grid))
    likelihood = array(0, dim = c(grid,grid))
    ll = array(0, dim = c(grid,grid))

    for (s in 1:length(sigma))
    {	
	for (m in 1:length(u))
	{	

		Px_g_u_s[[s,m]] = 1/(sqrt(2*pi)*sigma[s])*exp(-1*(gauss-u[m])^2/(2*(sigma[s]^2)))
		likelihood[s,m] = 1/(sqrt(2*pi)*sigma[s])*exp((-1*sum((data-u[m])^2))/(2*sigma[s]^2))
		plot(gauss,Px_g_u_s[[s,m]], type='l',axes=FALSE,ylim = c(0,4), frame.plot=F)
		axis (1,at = c(0,max(gauss)),labels=FALSE)
		
		#print(range(Px_g_u_s[[s,m]]))
	}
    }

     mtext("Gaussian Hypothesis Space", side = 3, line = -3,outer = TRUE)
   
    sigma = rev(sigma)
    for (s in 1:length(sigma))
    {	
        for (m in 1:length(u))
        {
            ll[s,m] = -n*log2(sqrt(2*pi)*sigma[s]) - (sum(((data - u[m])^2)/(2*sigma[s]^2)))
        }
    }
    dev.new()
    par(mfrow=c(1,1))
    contour(x=u,y=sigma,(likelihood),xlim=range(u),ylim = rev(range(sigma)),levels=pretty(range(likelihood),50),drawlabels=TRUE,axes=TRUE,main = "Log Likelihood",xlab="u", ylab="sigma")

    sigma=rev(sigma)


    dev.new()
    par(mfrow = c(grid,grid),mar =  c(.5,.5,.5,.5))


   
    max_likelihood = mean(data)

    r = range(likelihood)
    lwd = seq(from = r[1],to=r[2],by =(r[2]-r[1])/5)
    total=0
    for (s in 1:length(sigma))
    {	
	for (m in 1:length(u))
	{	
	
		posterior = Px_g_u_s[[s,m]]*likelihood[s,m]
		total =total+posterior
		if (likelihood[s,m] > (exp(-8)*max_likelihood))
		{
			w = which.min(abs(lwd -likelihood[s,m] ))
			plot(gauss,posterior, type='l',lwd=w,ylim = c(0,.2),xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
			axis (1,at = c(0,max(gauss)),labels=FALSE)
		}else
			plot.new()
	}
    }
    mtext("Gaussian Likelihood", side = 3, line = -3,outer = TRUE)
    
    total=total/(length(u)+length(sigma))
    return( total)

}



GaussianLikelihoodMixture <- function(data,u,sigma,mpi)
{
    gridu = 4
    grids = 5

    pi1 = mpi
    pi2 = 1-pi1
    
    sigmaf = sigma[1]
    sigmat = sigma[2]
    sigma1 = rev(seq(from = sigmaf, to = sigmat, by = (sigmat-sigmaf)/(grids-1)))
    sigma2 = rev(seq(from = sigmaf, to = sigmat, by = (sigmat-sigmaf)/(grids-1)))

    meanf = u[1]
    meant = u[2]
    u1 = seq(from = meanf, to = meant, by = (meant-meanf)/(gridu-1))
    u2 = seq(from = meanf, to = meant, by = (meant-meanf)/(gridu-1))

    gaussf = 0
    gausst = 3
    gauss = seq(from = gaussf, to = gausst, by = (gausst- gaussf)/90)

    Px_g_u_s = array(list(), dim = c(length(pi1),grids,grids,gridu,gridu))
    likelihood = array(0, dim = c(length(pi1),grids,grids,gridu,gridu))


  
    maxlikelihood = mean(data)

    for ( p in 1:length(pi1))
    {

	dev.new()
	par(mfrow = c(grids^2,gridu^2),mar = c(.5,.5,.5,.5))
          
	for(s1 in 1:length(sigma1))
	{

		for(s2 in 1:length(sigma2))
		{

		
			for (m1 in 1:length(u1))
			{


				for (m2 in 1:length(u2))
				{


	Px_g_u_s[[p,s1,s2,m1,m2]] = ((pi1[p]/(sqrt(2*pi)*sigma1[s1]))*exp(-1*(gauss-u1[m1])^2/(2*(sigma1[s1]^2)))) +

		((pi2[p]/(sqrt(2*pi)*sigma2[s2]))*exp(-1*(gauss-u2[m2])^2/(2*(sigma2[s2]^2))))

	likelihood[p,s1,s2,m1,m2] = (pi1[p]/(sqrt(2*pi*sigma1[s1]^2))*exp((-1*sum((data-u1[m1])^2))/(2*sigma1[s1]^2))) +

		(pi2[p]/sqrt((2*pi*sigma2[s2]^2))*exp((-1*sum((data-u2[m2])^2))/(2*sigma2[s2]^2)))

	plot(gauss,Px_g_u_s[[p,s1,s2,m1,m2]], type='l',xaxt='n',yaxt='n', ann=FALSE,  frame.plot=F)
	axis (1,at = c(0,max(gauss)),labels=FALSE)

	
				
				}
			}
		}

	}
 mtext(paste("Mixed Gaussian Hypothesis Space pi=",pi1[p]), side = 3, line = -3,outer = TRUE)
    	
    }
    
    total5=0
    r = range(likelihood)
    lwd = seq(from = r[1],to=r[2],by =(r[2]-r[1])/5)
    for ( p in 1:length(pi1))
    {

	dev.new()
	par(mfrow = c(grids^2,gridu^2),mar = c(.5,.5,.5,.5))

	for(s1 in 1:length(sigma1))
	{

		for(s2 in 1:length(sigma2))
		{

		
			for (m1 in 1:length(u1))
			{


				for (m2 in 1:length(u2))
				{

					posterior = Px_g_u_s[[p,s1,s2,m1,m2]] * likelihood[p,s1,s2,m1,m2]
					total5 = total5 +posterior
					if (likelihood[p,s1,s2,m1,m2] > (exp(-8)*maxlikelihood))
					{
						w = which.min(abs(lwd -likelihood[p,s1,s2,m1,m2] ))
						plot(gauss,posterior, type='l',xaxt='n',yaxt='n', ann=FALSE,lwd = w,  frame.plot=F)
						axis (1,at = c(0,max(gauss)),labels=FALSE)
					}else 
						plot.new()
				}
			}
		}
	}
        mtext(paste("Mixed Gaussian Likelihood pi =",pi1[p]), side = 3, line = -3,outer = TRUE)
   	
    }
   
    total5=total5/(length(pi1)*(grids^2+gridu^2))
    dev.new()
    par(mfrow=c(1,1))
    plot(total5, type='l',main = "Posterior Probability",ylab="")
    return (total5 )
}




data = c(.8,.75,.77,.78,1.8)

#range of u
u = c(0,2)

#range of sigma
sigma = c(.1, 1)
    
GaussianLikelihood(data,u,sigma)


mpi = c(.6)

GaussianLikelihoodMixture(data,u,sigma,mpi)
