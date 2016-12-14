PostProcess_beta <- function(MCMC_beta,model,wavespecsc,Ty,pcaspecsx,testspecs,wave_list)
{
  #POSTPROCESS_BETA	Performs post-processing in the x-wavelet space for the
  #                   function-on-function regression model described in
  #                   Meyer et al. 2014 (Biometrics). Requires inputs from
  #                   the WFMM exectuable availble from the Biostatistic
  #                   Group at MD Anderson.
  #
  #       INPUTS:
  #       MCMC_BETA   MCMC samples of beta with rows indicating sample
  #
  #       MODEL       The model structure resulting from calling wfmm_setup
  #
  #       WAVESPECSC  Wavelet specs for the x-space wavelets
  #       
  #       Ty          Grid length for y
  #
  #       PCASPECSX   PCA specs for the x-space
  #
  #       TESTSPECS   Specs for posterior functional inference
  #
  #       WAVE_LIST   Coefficients from DWT of each row of original X
  #
  #   Created: 9/9/2014
  #   By: Mark John Meyer
  
  ## Initialize output ##
  postout <- list()
  
  ## parameters ##
  p <- dim(model$X)[2]
  
  Tx <- wavespecsc$t
  delt <- testspecs$delt
  alf <- testspecs$alf
  twosample <- testspecs$twosample
  if(twosample == 1)
  {
    delt0 <- testspecs$delt[1]
    delt1 <- testspecs$delt[2]
  }
  
  ## unlist x-space PCA specs ##
  if(pcaspecsx$pca == 1)
  {
    pcalevelx <- pcaspecsx$pca_level;            
    scorex <- pcaspecsx$score;
    coefx <- pcaspecsx$coef;
    meansx <- pcaspecsx$mean_mat;
  }
  
  
  ##
  ##
  
  ## update x-space wavelet specs ##
  if(wavespecsc$compress == 1)
  {  wavespecsc$Kj <- wavespecsc$Kj_all}
  
  ## perform IDWT in X=space ##
  #  WFMM performs Y-space IDWT
  B <- dim(MCMC_beta)[1];
  if(twosample == 0)
  {
    Beta <- matrix(NaN,B,Tx*Tx)
    A <- matrix(NaN,B,Tx)
    for(i in 1:B)
    {
      ## separate b(t) and a(t) ##
      # Mbi <- MCMC_beta[i,] #-------------------------------------Check these lines
      # dim(Mbi) <- c(Ty,p) #------------------------------------------
      # abi <- t(Mbi) #-------------------------------------------------
       ## OR ##
      abi <- MCMC_beta[i,]
      dim(abi) <- c(Ty,p) #no transpose needed
      ai <- abi[1,]
      yidwt <- abi[2:dim(abi)[1],]
      
      ## update intercept ##
      A[i,] <- ai
      
      ## PCA and IDWT ##
      if(pcaspecsx$pca == 1)
      {
        pcaCol <- pcaspecsx$output[(pcaspecsx$output[,1] == pcalevelx),2]
        temp <- matrix(0,dim(scorex)[2],Tx)
        temp[1:pcaCol,] <- yidwt
        Betai <- t(t(temp)%*%solve(coefx) + meansx)          
      }
      if(wavespecsc$compress == 1)
      {
        temp <- matrix(0,dim(t(Betai))[1],length(wavespecsc$keep))
        temp[,wavespecsc$keep==1] <- t(Betai)
        Betai <- temp
      }
      Betai <- t(idwt_rows(Betai,wavespecsc))#------------------------------Uses function idwt_rows
      Beta[i,] <- as.vector(Betai) #--------------------------------------------------------------Double-check
      
      if(i%%100 == 0)
      {  cat('\n Done projecting MCMC samples M = ',i,'\n')}
    }
  }else{
    Beta0 <- matrix(NaN,B,Tx*Tx)
    Beta1 <- matrix(NaN,B,Tx*Tx)
    Diff <- matrix(NaN,B,Tx*Tx)
    A0 <- matrix(NaN,B,Ty)
    A1 <- matrix(NaN,B,Ty)
    for(i in 1:B)
    {
      #Betai = reshape(MCMC_beta(i,:)',Ty,p)';
      Mbi <- MCMC_beta[i,] #-------------------------------------Check these lines
      dim(Mbi) <- c(Ty,p) #------------------------------------------
      Betai <- t(Mbi) #-------------------------------------------------
                                                                                 
      ## split off surfaces ##
      As <- Betai[1:2,]
      cut_p <- dim(pcaspecsx$X)[2]
      yidwt0 <- Betai[3:(cut_p+2),]
      yidwt1 <- Betai[(cut_p+3):dim(Betai)[1],]
     
      ## remove and update intercepts ##
      A0[i,] <- As[1,]
      A1[i,] <- As[2,]
     
      ## inverse PCA ##
      if(pcaspecsx$pca==1)
      {
        ##  get # of PCA columns ##
        pcaCol <- pcaspecsx$output[pcaspecsx.output[,1]==pcalevelx,2]
     
        ## insert zeros for unused components ##
        temp0 <- matrix(0,dim(scorex)[2],Tx)
        temp1 <- matrix(0,dim(scorex)[2],Tx)
        temp0[1:pcaCol,] <- yidwt0
        temp1[1:pcaCol,] <- yidwt1;
     
        ## project back into data space ##
        beta_temp0 <- t(t(temp0)%*%solve(coefx) + meansx) #--------------------------------check dimensions here
        beta_temp1 <- t(t(temp1)%*%solve(coefx) + meansx) #-----------------------------------and here
      }else{
        beta_temp0 <- yidwt0
        beta_temp1 <- yidwt1
      }
      if(wavespecsc$compress== 1)
      {
        ## first surface ##
        temp0 <- matrix(0,dim(beta_temp0)[2],length(wavespecsc$keep))
        temp0[,wavespecsc$keep==1] <- t(beta_temp0)
        beta_temp0 <- temp0
        ## second surface ##
        temp1 <- matrix(0,dim(beta_temp1)[2],length(wavespecsc$keep))
        temp1[,wavespecsc$keep==1] <- t(beta_temp1)
        beta_temp1 <- temp1
      }
      ## x-space IDWT ##
      xidwt0 <- t(idwt_rows(beta_temp0,wavespecsc)) #-------------------------check dimensions from idw_rows
      xidwt1 <- t(idwt_rows(beta_temp1,wavespecsc)) #---------------------------------again
  
      ## update surfaces ##
      Beta0[i,] <- as.vector(xidwt0) #-------------------------Double-check this (dimensions,etc.)
      Beta1[i,] <- as.vector(xidwt1)
      Diff[i,]  <- Beta1[i,] - Beta0[i,]
  
      if((i%%100) == 0)
        cat('\n Done projecting MCMC samples M = ',i,'.\n')
    }
  }
  
  ## implement BFDR %%
  if(twosample == 0)
  {
    R <- dim(Beta)[2]
    pst <- rep(NaN,R)
    lFDR <- rep(NaN,R)
    for(r in 1:R)
    {
      Betar <- abs(Beta[,r])
      pst[r] <- (1/B)*length(which(Betar > delt))
      if(pst[r] == 1)
      {  
        pst[r] <- 1 -(2*B)^(-1)
      }
      lFDR[r] <- 1 - pst[r]
    }
    if(sum(pst) == 0)
    {  
      stop("No coef > delta")
    }
    pr <- sort(pst,decreasing = T)
    rstar <- cumsum(1-pr)/c(1:R)
    gam <- max(which(rstar <= alf))
    if(length(gam)==0)
    {
      phi <- 1
    }else{
      phi <- pr[gam]
    }
    psi <- (pst >= phi)
  }else{
    ## BFDR on Diff ##
    R  <- dim(Diff)[2]
    pst <- rep(NaN,R)
    lFDR <- rep(NaN,R)
    for(r in 1:R)
    {  
      Diffr <- abs(Diff[,r])
      pst[r] <- (1/B)*length(which(Diffr > delt))
      if(pst[r] == 1)
      {
        pst[r] <- 1 - (2*B)^(-1)
      }
      lFDR[r] <- 1 - pst[r]            
    }
    if(sum(pst) == 0)
    {
      stop("No coef > delta")
    }
    pr <- sort(pst,decreasing = T)
    rstar <- cumsum(1-pr)/c(1:R)
    gam <- max(which(rstar <= alf))
    if(length(gam)==0)
    {
      phi <- 1
    }else{
      phi <- pr[gam]
    }
    psi <- pst >= phi
    
    ## BFDR on Beta0 ##
    R0 <- dim(Beta0)[2]
    pst0 <- rep(NaN,R0)
    lFDR0 <- rep(NaN,R0)
    for(r in 1:R0)
    {
      Betar0 <- abs(Beta0[,r])
      pst0[r] <- (1/B)*length(which(Betar0 > delt0))
      if(pst0[r] == 1)
      {
        pst0[r] <- 1 - (2*B)^(-1)
      }
      lFDR0[r] <- 1 - pst0[r]
    }
    if(sum(pst0) == 0)
    {
      stop("No coef > delta")
    }
    pr0 <- sort(pst0,decreasing = T)
    rstar0 <- cumsum(1-pr0)/c(1:R0)
    gam0 <- max(which(rstar0 <= alf/2))
    if(length(gam0)==0)
    {
      phi0 <- 1
    }else{
      phi0 <- pr0[gam0]
    }
    psi0 <- (pst0 >= phi0)
    
    ## BFDR on Beta1 ##
    R1 <- dim(Beta1)[2]
    pst1 <- rep(NaN,R1)
    lFDR1 <- rep(NaN,R1)
    for(r in 1:R1)
    {  
      Betar1 <- abs(Beta1[,r])
      pst1(r) <- (1/B)*length(which(Betar1 > delt1))
      if(pst1[r] == 1)
      {
        pst1[r] <- 1 - (2*B)^(-1)
      }
      lFDR1(r) <- 1 - pst1[r]
    }
    if(sum(pst1) == 0)
    {
      stop("No coef > delta")
    }
    pr1 <- sort(pst1,decreasing = T)
    rstar1 <- cumsum(1-pr1)/c(1:R1)
    gam1 <- max(which(rstar1 <= alf/2))
    if(length(gam1)==0)
    {
      phi1 <- 1
    }else{
      phi1 <- pr1[gam1]
    }
    psi1 <- (pst1 >= phi1)
  }
      
  dim(psi) <- c(Tx,Tx)
  postout$psi <- psi
  cat("\n Done calculating FDR.\n \n")
      
  ## average samples and find credible interval, update postout ##
  if(twosample == 0)
  {
    Q025_bhat <- apply(Beta,2,quantile,probs = 0.025, type = 5) #---------------------------------------Quantiles calculated slightly differently? Look up method compared to Matlab
    Q975_bhat <- apply(Beta,2,quantile,probs = 0.975, type = 5) #-----------------------------------------------again
    pwCI <- rep(1,length(Q975_bhat))
    for(i in 1:length(pwCI))
    {
      if(Q025_bhat[i] < 0 & Q975_bhat[i] > 0)
      {
        pwCI[i] <- 0
      }
    }
  
    ## beta credible intervals and estimate ##
    dim(pwCI) <- c(Tx,Tx); dim(Q025_bhat) <- c(Tx,Tx); dim(Q975_bhat) <- c(Tx,Tx)
    meanBeta <- colMeans(Beta); dim(meanBeta) <- c(Tx,Tx)
    postout$pwCI <- pwCI
    postout$Q025_bhat <- Q025_bhat
    postout$Q975_bhat <- Q975_bhat
    postout$bhat <- meanBeta
    
    ## intercept credible intervals and estimate ##
    postout$Q025_ahat <- apply(A,2,quantile,probs = 0.025,type = 5)
    postout$Q975_ahat <- apply(A,2,quantile,probs = 0.975,type = 5)
    postout$ahat <- colMeans(A)
  }else{
    ## group 0 ##
    Q025_bhat0 <- apply(Beta0,2,quantile,probs = 0.025,type = 5)
    Q975_bhat0 <- apply(Beta0,2,quantile,probs = 0.975,type = 5)
    pwCI0 <- rep(1,length(Q975_bhat0))
    for(i in 1:length(pwCI0))
    {
      if(Q025_bhat0[i] < 0 & Q975_bhat0[i] > 0)
      {
        pwCI0[i] <- 0
      }
    }
  
    ## beta0 credible intervals and estimate ##
    dim(pwCI0) <- c(Tx,Tx); dim(Q025_bhat0) <- c(Tx,Tx); dim(Q975_bhat0) <- c(Tx,Tx)
    meanBeta0 <- colMeans(Beta0); dim(meanBeta0) <- c(Tx,Tx)
    postout$pwCI0 <- pwCI0
    postout$Q025_bhat0 <- Q025_bhat0
    postout$Q975_bhat0 <- Q975_bhat0
    postout$bhat0 <- meanBeta0
  
    ## group 1 ##
    Q025_bhat1 <- apply(Beta1,2,quantile,probs = 0.025,type = 5)
    Q975_bhat1 <- apply(Beta1,2,quantile,probs = 0.975,type = 5)
    pwCI1 <- rep(1,length(Q975_bhat1))
    for(i in 1:length(pwCI1))
    {  
      if(Q025_bhat1[i] < 0 & Q975_bhat1[i] > 0)
      {
        pwCI1[i] <- 0
      }
    }
  
    ## beta1 credible intervals and estimate %%
    dim(pwCI1) <- c(Tx,Tx); dim(Q025_bhat1) <- c(Tx,Tx); dim(Q975_bhat1) <- c(Tx,Tx)
    meanBeta1 <- colMeans(Beta1); dim(meanBeta1) <- c(Tx,Tx)
    postout$pwCI1 <- pwCI1
    postout$Q025_bhat1 <- Q025_bhat1
    postout$Q975_bhat1 <- Q975_bhat1
    postout$bhat1 <- meanBeta1
  
    ## diff ##
    Q025_diff <- apply(Diff,2,quantile,probs = 0.025,type = 5)
    Q975_diff <- apply(Diff,2,quantile,probs = 0.975,type = 5)
    pwCId <- rep(1,length(Q975_diff)) #-----------------------------Q975_bhatd doesn't exist???? Double-check this
    for(i in 1:length(pwCId))
    {  
      if(Q025_diff[i] < 0 & Q975_diff[i] > 0)
      {
        pwCId[i] <- 0
      }
    }
    
    ## diff credible intervals and estimate ##
    dim(pwCId) <- c(Tx,Tx); dim(Q025_diff) <- c(Tx,Tx); dim(Q975_diff) <- c(Tx,Tx)
    meanDiff <- colMeans(Diff); dim(meanDiff) <- c(Tx,Tx)
    postout$pwCId <- pwCId
    postout$Q025_diff <- Q025_diff
    postout$Q975_diff <- Q975_diff
    postout$diff <- meanDiff
    
    ## intercepts credible intervals and estimates ##
    postout$Q025_ahat0 <- apply(A0,2,quantile,probs = 0.025,type = 5)
    postout$Q975_ahat0 <- apply(A0,2,quantile,probs = 0.975,type = 5)
    postout$ahat0 <- colMeans(A0)
    
    postout$Q025_ahat1 <- apply(A1,2,quantile,probs = 0.025,type = 5)
    postout$Q975_ahat1 <- apply(A1,2,quantile,probs = 0.975,type = 5)
    postout$ahat1 <- colMeans(A1)
    
    ## group specific psi ##
    dim(psi0) <- c(Tx,Tx); dim(psi1) <- c(Tx,Tx)
    postout$psi0 <- psi0
    postout$psi1 <- psi1
  }
  cat("\n Done averaging and finding CIs.\n \n")
  
  ## calculate SimBaS ##
  if(twosample == 0)
  {
    bandmaps <- jointband_maps(Beta,alf) #---------------------------------------Check jointband_maps function output
    dim(bandmaps$SimBaS) <- c(Tx,Tx); dim(bandmaps$upper_CI) <- c(Tx,Tx); dim(bandmaps$lower_CI) <- c(Tx,Tx)
    postout$SimBaS <- bandmaps$SimBaS
    postout$USimBaS <- bandmaps$upper_CI
    postout$LSimBaS <- bandmaps$lower_CI
  }else{
    ## group 0 ##
    bandmaps0 <- jointband_maps(Beta0,alf/2) #---------------------------------------------Check jointband_maps output
    dim(bandmaps0$SimBaS0) <- c(Tx,Tx); dim(SimBaS0$upper_CI) <- c(Tx,Tx); dim(SimBaS0$lower_CI) <- c(Tx,Tx)
    postout$SimBaS0 <- bandmaps0$SimBaS0
    postout$USimBaS0 <- bandmaps0$upper_CI0
    postout$LSimBaS0 <- bnadmaps0$lower_CI0
    
    ## group 1 ##
    bandmaps1 <- jointband_maps(Beta1,alf/2)
    dim(bandmaps1$SimBaS1) <- c(Tx,Tx); dim(SimBaS1$upper_CI) <- c(Tx,Tx); dim(SimBaS1$lower_CI) <- c(Tx,Tx)
    postout$SimBaS1 <- bandmaps1$SimBaS1
    postout$USimBaS1 <- bandmaps1$upper_CI1
    postout$LSimBaS1 <- bandmaps1$lower_CI1
    
    ## diff ##
    bandmapsD <- jointband_maps(Diff,alf)
    dim(bandmapsD$SimBaSD) <- c(Tx,Tx); dim(SimBaSD$upper_CI) <- c(Tx,Tx); dim(SimBaSD$lower_CI) <- c(Tx,Tx)
    postout$SimBaSD <- bandmapsD$SimBaSd
    postout$USimBaSD <- bandmapsD$upper_CId
    postout$LSimBaSD <- bandmapsD$lower_CId
  }
  
  fprintf('\n Done calculating MAPs.\n \n');
  return(list("postout"=postout))
}