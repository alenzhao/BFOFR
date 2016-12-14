#' Discrete Wavelet Transformation and Compression
#'
#' @export
#' @description This function performs a DWT on each row (observed function) of X according to methods described in Morris and Carroll (2006).
#' In addition, performs joint compression on all of the rows of wavelet coefficients as done in Morris et al. (2011). Note that whereas in Meyer et al. (2015)
#' this function can accomdate multi-dimensional wavelet decomposition, our function only works for 1-dimension.
#' @param X_raw The original observed matrix X of functions sample over a grid.
#' @param wavespecs A list containing the specifications for the DWT to be performed on X_raw.
#' @param alpha compression levels; i.e, P=0.95 represents a compression ratios of over 100:1, which is considered as compression levels. See more details in Morris et al,(2011).
#' @param xticks compression paraemters for X. See more details in Morris et al, (2011).
#' @return Output from this function is a list containing 1) "compRes", 2) "compSet", 3)"D_all", 4)"C", 5) "keep", and 6) "wavespecsc".
#' CompRes: return compression result; CompSet: settings for compression on X; D: wavelet coefficients for X;  C: wavelet levels; Keep: an indicator vector with 1s and 0s.
#' 1s mean we retain its correspoinding columns in the compression and 0s mean we exclude its correspoinding columns after the compression.
#' @examples
#'

Wavelet_compress<- function(X_raw,wavespecs,alpha,xticks)
{
  ## Dependencies
  library(wavelets,quietly = T)
  #library(wavethresh,quietly = T)

  tick <- proc.time()[3]

  # 1203 comment the below conditions out since those conditions are dependent on the # of inputs, which we have 6 inputs.
  # If it doesn't work, we can retrieve it.
  # if(nargs()<3)
  #   Image_matrix <- 1
  #
  # if(nargs()<2)
  # {
  #   wavespecs <- list()
  #   wavespecs$wavelet<-'db4';
  #   wavespecs$wtmode<-'symmetric';
  #   wavespecs$nlevels<-floor(log(min(dim(X_raw)))/log(2)); #check this result compared to Matlab
  #   wavespecs$ndim<-2;
  # }
  #
  # if(nargs()<4)
  #   alpha<-c(.9,.95,.99,.995,.999)
  #
  # if(nargs()<5)
  #   xticks<-c(200,300,400,500,1000,2000,3000,4000,5000,7500,10000,20000)
  #
  # if(nargs()<6)
  #   graph <- 0

  x_min<-min(alpha);
  x_max<-.999;
  x_step<-0.001;
  x_length<-floor((x_max-x_min)/x_step)+1;
  x_tick<-(0:x_length)*x_step+x_min;

  x_tick<-c(.95+(0:49)*.001,(.999+(1:9)*.0001)); #-------------keep this or the other x_tick? x_tick is fine.
  yticks<-c(0.95,0.975,0.99,0.995,0.999,0.9999);
  x_length<-length(x_tick)-1; #--------------------------------keep this or the other x_length? x_length is fine.
  # We are dealing with 1-D dimensional case.And we will make our codes more flexible, which can be used to high dimensional case.
  if(wavespecs$ndim==1){
    n <- dim(X_raw)[1]
    X <- X_raw[1,]

    # Initalize list in which to store wavelet transforms of each row
    wave_list <- list()
    names_list <- rep(NA,n)

    # 1203 help(dwt) check this "boundary = wavespecs$wtmode". We change it to "periodic"
    wave <- dwt(X,filter = wavespecs$wavelet,n.levels = wavespecs$nlevels,boundary ="periodic") #need to change filter and boundary specs to fit with the R function###############check the R code
    C <- c(unlist(wave@W),unlist(wave@V[length(wave@V)])) ## 224 not 225 12/07
    D <- matrix(0,n,length(C));
    S <- rep(0,wavespecs$nlevels+2)
    S[1] <- length(wave@V[[length(wave@V)]])
    for(i in 2:(wavespecs$nlevels+1))
    {
      S[i] <- length(wave@W[[i-1]])
    }
    S[length(S)] <- sum(S[1:length(S)-1])
    wavespecs$Kj <- S[1:(length(S)-1)];
    wavespecs$t <- S[length(S)];#################### 12/07
    for(i in 1:n)
    { #1203 same issue with setup for boundary. It was "boundary = wavespecs$wtmode". Need to check it again.
      # Error in wt.filter(filter) : Invalid filter name.
      wt<-dwt(X_raw[i,],filter = wavespecs$wavelet,n.levels = wavespecs$nlevels,boundary = "periodic")
      wave_coef <- list(W=wt@W,V=wt@V)
      wave_list[[i]] <- wave_coef
      names_list[i] <- paste("wave",as.character(i),sep = "")
      # newly added
      D[i,] <-c(unlist(wt@W),unlist(wt@V[length(wt@V)]))
      #S[1] <- length(wt@V[[length(wt@V)]])
      #for(i in 2:(wavespecs$nlevels+1))
      #{
      #  S[i] <- length(wt@W[[i-1]])
      #}
      # newly added
    }
    names(wave_list) <- names_list
    wave_list$nlevels <- wt@level
    wave_list$boundary <- wt@boundary
    wave_list$filter <- wavespecs$wavelet
  }
  cat('Done with Wavelet Transforms\n')
  tock <- proc.time()[3]-tick
  tock
  #############################################1127###########
  #rm(S,X,X_raw)

  En <- D #single(D);
  total_Energy <- matrix(0,1,dim(D)[1]);
  for (i in 1:n){
    Csort <-sort(-abs(D[i,]),index.return=T);
    total_Energy[i] <- sum(Csort$x^2);''
    energy <-cumsum(Csort$x^2)/total_Energy[i]
    En[i,Csort$ix] <- energy
    i;
  }
  cat('Done Computing Relative Energy\n')


  Ncoeffs<-matrix(0,dim(D)[1],x_length+1);
  for(j in 1:(x_length+1)){
    temp<-sum(En<x_tick[j])
    for(i in 1:n){
      Ncoeffs[i,j]<-sum(temp>(i-1));
    }
    j;
  }

  cat('Done Computing # coefficients,', proc.time()[3]-tick, "seconds\n")

  En_min<-matrix(0,dim(Ncoeffs)[1],dim(Ncoeffs)[2]);
  for(j in 1:(x_length+1))
  {
    temp<-colSums(En<x_tick[j])
    for(i in 1:n)
    {
      keep<-(temp>(i-1))
      if(sum(keep)>1)
      {
        En_min[i,j]<-min(colSums(t(D[,keep])^2)/total_Energy)
      }
      if(sum(keep)==1)
      {
        En_min[i,j]<-min(t(D[,keep])^2/total_Energy)
      }
      if(sum(keep)==0)
      {
        En_min[i,j]<-0;
      }
    }
    j;
  }

  cat('Done Computing Minimum Energies', proc.time()[3]-tick, "seconds\n")

  ##### Jonathan picks up from here
  # testing with Matlab codes line by line
  # library(R.matlab,quietly = T)
  #setwd('C:/Users/Zhaohu/Dropbox/2016FDA/TestingDataset')
  #library(R.matlab)
  #x_tick<- readMat('x_tick.mat') # input x_tick pass
  #names(x_tick)=c()
  #x_tick<-unlist(x_tick)
  #x_tick<-matrix(x_tick, nrow=1, byrow=TRUE)
  settings<-cbind(x_tick,x_tick,x_tick,x_tick)  # setting pass # 59*4
  rows<-1:n;
  #need to introduce an new variable ; remove it when it is necessary.
  for(k in 1:length(x_tick)){
    min_energy<-x_tick[k];

    (t(En_min)>min_energy)
    ff1<-ifelse((t(En_min)>min_energy),1,0)
    ff2<-colSums(ff1, na.rm=FALSE)
    ff3<-ifelse(ff2>0,1,0)
    i<-which(rows[ff3]>0, arr.ind = TRUE)
    i2 <- rows[1:dim(En_min)[1]][colSums(t(En_min)>min_energy)>0]

    hu1<-t(En_min[i,])>min_energy
    hu2<-ifelse((t(En_min[i,])>min_energy),1,0)
    hu3<-colSums(hu2, na.rm=FALSE)
    j<-length(x_tick)-hu3
    j2 <- (length(x_tick)+1)-colSums(t(En_min[i2,])>min_energy)


    #Ncoeffs<- readMat('Ncoeffs.mat') # input Ncoeffs pass
    #names(Ncoeffs)=c()
    #Ncoeffs<-unlist(Ncoeffs)
    #Ncoeffs<-matrix(Ncoeffs, nrow=200, byrow=TRUE)

    #compression<-(diag(Ncoeffs(i,j)));
    compression<-(diag(Ncoeffs[i2,j2]));
    #istar<-min(rows(compression==min(compression))); find the index of the first smallest value of compression.
    istar<- which(compression==min(compression))
    istar<-min(istar)


    settings[k,3]<-istar;
    #cat("Settings 3\n")
    settings[k,4]<-x_tick[j2[istar]];  #1204 Issue here
    #cat("Settings 4\n")
    settings[k,2]<-compression[istar];# pass
    #cat("settings 2\n")
    #cat(k,"\n")
  }

  #cat("Done creating settings")

  ## 12-3-make visulization
  ncoeffs<-settings[,2];
  # if(graph==1){
  #   plot(ncoeffs,settings[,1],main = "Minimum % Energy Preserved vs. Number of Coefficient",xlab="Number of Coefficient",ylab = "Minimum  %Energy Preserved")
  # }
  alpha_row<-matrix(0,length(alpha),1);
  for (i in 1:length(alpha)) {
    alpha_row[i]<-sum(settings[,1]<alpha[i])+1;
  }
  #cat("Done creating alpha_row matrix")

  proc.time()[3]-tick

  #En<- readMat('En.mat') # input En pass
  #names(En)=c()
  #En<-unlist(En)
  #En<-matrix(En, nrow=200, byrow=TRUE)


  #D<- readMat('D.mat') # input En pass
  #names(D)=c()
  #D<-unlist(D)
  #D<-matrix(D, nrow=200, byrow=TRUE)

  results<-settings[alpha_row,];
  keep<-matrix(0,dim(results)[1],dim(D)[2]);
  for(i in 1:dim(results)[1]){
    # i=1; results[i,4];results[i,4]
    # keep[i,]<-sum(En<results[i,4])>=results[i,3]==1; #-----------------------CHECK THIS LINE; checked
    # match with the Matlab code
    zz1<-ifelse(En<results[i,4],1,0)
    zz1<-colSums(zz1, na.rm=FALSE)
    keep[i,]<-ifelse(zz1>results[i,3],1,0)
    #match with the Matlab code
  }


  # function output

  waveletresult<-list("compRes"=results,"compSet"=settings,"D_all"=D, "C"=C,"keep"=keep,"wavespecsc"=wavespecs)
  return(waveletresult)
}

