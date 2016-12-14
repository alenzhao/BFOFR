#' Function to generate beta surfaces
#'
#' @export
#' @description This function generates one of 3 beta surfaces examined in the simulations of Meyer, et. al. (2015): 1) ridge, 2) lagged, and 3) immediate.
#' The bivariate function is evaluated over grids t and v, assumed to be of equal length and with equal spacing. The length of the grid must be
#' specified. A 3D perspective plot of the generated surface may also be output.
#' @param surface Choose either "ridge", "lagged", or "immediate".
#' @param len The length of the grids over which the surface is evaluated.
#' @param surfplot An indicator of whether or not to output a 3D-perspective plot of the surface. The default is False.
#' @return Output is a matrix of values corresponding to the bivariate function evaluated over the gridpoints of t and v.
#' @examples
#' beta_gen("ridge",len = 225,F)

beta_gen=function(surface = "ridge", len, surfplot = F){
  if(!is.numeric(len))
  {
    stop("len must be a positive numeric (integer) value.\n")
  }

  ## Simulated Beta surfaces

  # grid length
  T. <-len

  # grids for simulation surfaces
  v <- matrix(data = seq(0,(T.-1),by = 1),nrow = T.,ncol = T.,byrow = T)
  t <- t(matrix(data = seq(0,(T.-1),by = 1),nrow = T.,ncol = T.,byrow = T))

  ## beta coefficients:

  if(surface=="ridge"){
    ## 1) ridge
    stnrs <- 1.75 # scale to control stnr
    b1var <- 0.003
    b1 <- stnrs*(1/(125*sqrt(2*pi*b1var)))*exp(-1/(2*b1var)*(t/T.-v/T.)^2)
  }else if(surface=="lagged"){
    ## 2) lagged
    b1var <- 0.003
    b1 <- (437/10000)*(1/sqrt(2*pi*b1var))*exp(-1/(2*b1var)*(t/T.-v/T.-0.5)^2)
  }else if(surface=="immediate"){
    ## 3) immediate
    b1 <- (1029/10000)*(1+1/(1+exp((0.25-t/T.+v/T.)/0.05)))
  }#else if(surface=="peak"){
  #  ## 4) peak
  #  sigma <- matrix(c(15, 0.5*15^2, 0.5*15^2, 15),2,2,byrow = T)
  #  a <- c(150,60)
  #  b1 <- (1225/10)*(1/(2*pi))*(det(sigma))^(-1/2)*exp(-t(a)%*%solve(sigma)%*%a/2)
  #}
  else{
    stop("Incorrect surface option.\n")
  }

  if(surfplot){
  # plot ridge
  zmin <- min(b1)
  zmax <- max(b1)
  persp(x = c(0:(T.-1)),y = c(0:(T.-1)),z = b1,zlim = c(zmin,zmax),col = "red",xlab = "X",ylab = "Y",zlab = "b1")
  #library(plotly,quietly = T)
  #plot_ly(x = c(0:224),y=c(0:224),z=b1) %>% add_surface()
  }

  return(b1)
}
