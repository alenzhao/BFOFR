#' Function to setup the design matrix
#'
#' @export
#' @description This function performs the specified basis operations on X (e.g. discrete wavelet transformation, PCA). This function should be used prior to using the wfmm executable file, which performs the BFOFR as described in Meyer, et. al. (2015)
#' @param Y N-by-T matrix of observed response functions. N is the total number of observed curves, T is the number of grid points at which the curves were observed.
#' @param X N-by-V matrix of observed predictor functions. N is the total number of observed curves, V is the number of grid points at which the curves were observed.
#' @param Z N-by-n random effects design matrix. N is the total number of observed curves, n is the number of subjects.
#' @param MCMCspecs A list of specifications for the MCMC algorithm used to estimate the beta coefficients. Needed for the wfmm executable file.
#' @param wavespecsx A list of specifications for the wavelet transformation of X. This list must include: "comp.level", a proportion supplied as the compression level, "alpha", "xticks", "nlevels", "wavelet", "wtmode", and "ndim"
#' @param pca An indicator of whether or not PCA is to be performed on X after the DWT. 1 indicates yes, 0 no. Default is 1.
#' @param testspecs A list of test specifications for posterior inference. This list must contain "twosample", an indicator of whether the design matrix has two samples, and FDR, a vector containing cutoffs for the Bayesian FDR procedure.
#' @param pcaspecsx If pca = 1, a list of specifications for the PCA performed on X. The list should contain "pca_level", the percentage of explained variation chosen as a cutoff to determine
#' how many PC's to keep, and "pca_cutoff", a vector of such pca level values to consider. If the argument is not supplied, defaults to
#' "pca_cutoff" = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9999) and "pca_level" = 0.99.
#' @return This function returns a list of 3 items: "model", "xspecs", and "testspecs". "model" contains "Z" (unchanged), "X" which is now the basis-space design matrix, and "C". "xspecs" is a list including the wavespecs and pcaspecs used for x. "testspecs" includes the original "testspecs" as well as splits the "FDR" vector into "alf" and "delt".
#' @examples
#'

wfmm_setup<-function(Y,X,Z,MCMCspecs,wavespecsx,pca,testspecs,pcaspecsx)
{
  ############################################################################
  # WFMM_SETUP Builds the design matrix for performing function-on-function  #
  #           regression as in Meyer et al. 2014 (Biometrics) using the wfmm #
  #           exectuable file availble from the Biostatistic Group at MD     #
  #           Anderson                                                       #
  #                                                                          #
  #                                                                          #
  #   Created: 7/18/2014                                                     #
  #   By: Mark John Meyer                                                    #
  ############################################################################

  # assume one sample if not otherwise specified
  if(nargs()<7)
    testspecs$twosample <- 0

  # default to PCA-based approach
  if(nargs()<6)
    pca <- 1

  # check and assign two-sample specs
  twosample <- testspecs$twosample
  if(twosample==1 && length(names(testspecs))<2)
    stop("Grouping variable required for two-sample model.")
  if(twosample==1)
    group <- testspecs$group

  # center and scale X
  N <- dim(X)[1]
  t <- dim(X)[2]
  scaleX <- scale(x = X,center = T,scale = T)

  # check if random effects exist, set sampleU accordingly
  if(all(is.na(Z))){
    MCMCspecs$SampleU <- 0
  }else{
    MCMCspecs$SampleU <- 1
  }

  # extract alpha levels and xticks for compression
  alpha <- wavespecsx$alpha
  xticks <- wavespecsx$xticks

  # compression performed as in Morris et al. ()
  wc <- Wavelet_compress(scaleX, wavespecsx, alpha, xticks)

  # select compression level
  comp_level <- wavespecsx$comp.level

  # additional compression output
  compOut <- list("compRes" = wc$compRes, "compSet" = wc$compSet, "C" = wc$C)

  # update wavespecs
  wavespecsc <- wc$wavespecsc #####
  wavespecsc$compOut <- compOut
  wavespecsc$keep <- wc$keep[alpha==comp_level,]
  wavespecsc$wave_list <- wc$wave_list
  wave_keep <- which(wc$keep[alpha==comp_level,]==1)
  wavespecsc$compress <- 1
  D <- wc$D_all[,wavespecsc$keep==1]

  #source("get_Kj_compress.R")
  # get new Kj
  g <- get_Kj_compress(wavespecsc)
  wavespecsc$Kj_all <- g$Kj_all
  wavespecsc$J_all <- length(g$Kj_all)
  wavespecsc$J <- length(g$Kj_comp)
  wavespecsc$K <- sum(g$Kj_comp)

  # perform PCA on X after DWT
  if(pca==1)
  {
    # Default PCA specs for X-specs
    if(nargs()<8)
    {
      pcaspecsx <- list("pca_cutoff" = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9999),"pca_level" = 0.99)
    }
    pca_cutoff = pcaspecsx$pca_cutoff
    princomp_D <- prcomp(D,center = T)
    coefs <- princomp_D$rotation
    score <- princomp_D$x
    eigs <- (princomp_D$sdev)^2 # technically these are only proportional to the eigenvalues
    comp_var <- cumsum(eigs)/sum(eigs)
    pca_keep <- rep(0,length(pca_cutoff))
    for(i in 1:length(pca_cutoff))
    {
      pca_keep[i] <- min(which(comp_var > pca_cutoff[i]))
    }

    pca_out <- cbind(pca_cutoff,pca_keep)

    # select PCA level
    pca_level <- pcaspecsx$pca_level

    # select PCs to use
    X <- score[,1:pca_out[(pca_out[,1] == pca_level),2]]

    # generate column mean matrix
    col_means <- mean(D)
    col_mean_mat <- matrix(0,nrow = t,ncol = dim(D)[2])
    for(i in 1:t)
    {
      col_mean_mat[i,] <- col_means
    }


    # save PCA specs
    pcaspecsx$X <- X;
    pcaspecsx$pca <- 1;
    pcaspecsx$score <- score;
    pcaspecsx$coef <- coefs;
    pcaspecsx$eigen <- eigs;
    pcaspecsx$mean_mat <- col_mean_mat;
    pcaspecsx$output <- pca_out;
    pcaspecsx$keep <- pca_out[(pca_out[,1] == pca_level),2]
  }else{
    pcaspecsx$pca <- 0
  }


  # twosample design matrix
  if(twosample == 1)
  {
    if(pca == 1)
    {
      # set up big design matrix
      X0 <- X[group == 0,]
      X1 <- X[group == 1,]
      twoSize <- dim(X0)[2] + dim(X1)[2]
      X <- matrix(0,nrow = N,ncol = twoSize)

      # fill in big design matrix
      X[group == 0, 1:dim(X0)[2]] <- X0;
      X[group == 1, (dim(X0)[2]+1):(twoSize)] <- X1;
    }else{
      # set up big design matrix
      tsX <- matrix(0,nrow = N,ncol = 2*dim(D)[2]);
      group0 <- which(group == 0);
      group1 <- which(group == 1);

      # fill in big design matrix
      tsX[group0, 1:dim(D)[2]] <- D[group0,]
      tsX[group1, (dim(D)[2]+1):(2*dim(D)[2])] <- D[group1,]

      # assign tsX to D
      D <- tsX
    }
  }

  # model Set up
  model <- list()
  if(twosample == 1)
  {
    # allows for group-specific intercepts
    a0 <- rep(0,N)
    a1 <- rep(0,N)
    a0[group == 0] <- 1
    a1(group == 1) <- 1

    if(gbf == 1) #----------------------------------------what is gbf???-----------------------------------------------------------------
    {
      model$X <- cbind(a0,a1,X)
    }else{
      model$X <- cbind(a0,a1,D)
    }
  }else{
    # include intercept
    a <- rep(1,N)
    if(pca == 1)
    {
      model$X <- cbind(a,X)
    }else{
      model$X <- cbind(a,D)
    }
  }
  model$Z <- Z
  model$C <- rep(1,dim(Y)[1])
  # model.m = 1; edit for more levels of random effects
  # model.H = 1;
  # model.Hstar = 0;

  # set up FDR for post-processing
  FDR <- testspecs$FDR;
  if(length(FDR) == 1)
  {
    testspecs$delt <- FDR[1];
    testspecs$alf <- 0.05;
  }
  if(length(FDR) == 2)
  {
    testspecs$delt <- FDR[1];
    testspecs$alf   = FDR[2];
  }
  if(length(FDR) == 3)
  {
    testspecs$delt <- FDR[1:2];
    testspecs$alf   = FDR[3];
  }
  if(length(FDR) > 3)
  {
    stop('Too many FDR inputs')
  }

  # function output
  xspecs <- list()
  xspecs$wavespecsc <- wavespecsc;
  xspecs$pcaspecs <- pcaspecsx;
  xspecs$wave_keep <- wave_keep

  return(list("model"=model,"xspecs"=xspecs,"testspecs"=testspecs))
}

