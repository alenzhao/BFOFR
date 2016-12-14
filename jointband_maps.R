jointband_maps <-function(MCMC_P,alpha)
{
  # function to calculate (1-alpha) joint credible band and 
  # multiplicity-adjusted probability score (MAPS) based on MCMC samples 
  # contained in a B*T matrix 'MCMC_P'. For details on the joint credible
  # band, refer to Crainiceanu et al 2007. At each location, 
  # the MAPS is defined as the minimum significance level at which the 
  # joint credible bands excludes zero. 
  
  # Inputs
  #    'MCMC_P' - a B by T matrix containing MCMC samples. 
  #        B - number of MCMC iterations (MCMCspecs.B)
  #        T - number of function samples (number of columns in Y)
  #    'alpha' - a row vector containing significance levels at which 
  #        joint credible bands are to be returned. 
  
  # Outputs
  #   'MAPS' - multiplicity-adjusted probability scores. MAPS are
  #    truncated at either 0.5 or 0.001 for reducing computational burden.
  #   'upper_CI' - a length(alpha) by T matrix containing the upper bounds 
  #        of the joint credible bands. The first row corresponds to 
  #        the first level in alpha.       
  #   'lower_CI' - a length(alpha) by T matrix containing the lower bounds 
  #        of the joint credible bands. 
  
  B <- dim(MCMC_P)[1]
  t <- dim(MCMC_P)[2]
  sd_P = rep(NaN,t)
  for(i in 1:t)
  {
    sd_P[i] <- sd(MCMC_P[,i])
  }
  mean_P <- colMeans(MCMC_P)
  z_P <- rep(NaN,B)
  for(j in 1:B)
  {
    z_P[j] <- max(abs((MCMC_P[j,]-mean_P)/sd_P))
  }
  levels <- sort(c(seq(0.01,0.5,0.01),seq(0.001,0.009,0.001),seq(0.0001,0.0009,0.0001)),decreasing = T)
  cb1 <- quantile(z_P,probs = (1-levels),type = 5)
  MAPS <- rep(0.5,t)
  for(k in 2:length(levels))
  {
    temp_ind <- (mean_P - cb1[k]*sd_P)*(mean_P + cb1[k]*sd_P)
    MAPS[temp_ind>0] <- levels[k]
  }
  
  if(nargs()>1)
  {
    m <- length(alpha)
    upper_CI <- matrix(0,m,t)
    lower_CI <- matrix(0,m,t)
    cb2 <- quantile(z_P,probs = (1-alpha),type = 5)
    for(j in 1:m)
    {
      upper_CI[j,] <- mean_P + cb2[j]*sd_P
      lower_CI[j,] <- mean_P - cb2[j]*sd_P
    }
  }
  return(list("MAPS"=MAPS,"upper_CI"=upper_CI,"lower_CI"=lower_CI))
}