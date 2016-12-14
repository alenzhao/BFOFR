#' Function to get Kj from compressed X
#'
#' @export
#' @description Given wavespecsc object with Kj = # coeffs per level and keep = vector of 0's and 1's with 1 in coordinate for coefficients to "keep",
#' figures out J  = # levels and Kj = # coefficients per level for the reduced set of coefficients that are kept after compression.
#' @param wavespecsc List of specifications for the wavelet transform of X. Specifically, wavespecsc$Kj, the number of coefs per level, and
#' wavespecsc$keep, are needed.
#' @return This outputs a list of two values: "Kj_comp" which is # of coefficients per level after the compression, and "Kj_all" which is
#' # of coefficients per level after the original wavelet transformation.
#' @examples
#'

get_Kj_compress <- function(wavespecsc)
{
  # Given wavespecs object with Kj=# coeffs per level and
  # keep = vector of 0's and 1's with 1 in coefficients to "keep"
  # figure out J and Kj = # levels and # coefficients per level for reduce
  # set of coefficients that are kept

  J=length(wavespecsc$Kj)
  Kj_all=wavespecsc$Kj
  Kj=rep(0,J)
  temp_first=c(1,1+cumsum(Kj_all[1:(length(Kj_all)-1)]))
  temp_last=cumsum(Kj_all)
  for (i in 1:J)
  {
    Kj[i]=sum(wavespecsc$keep[temp_first[i]:temp_last[i]])
  }
  # Also need to add what to do when there are "0 levels" or levels with
  # very small number of wavelet coefficients which we want to potentially
  # combine together
  #
  # Also, need to think about implications on spike-slab prior since we are
  # eliminating a lot of the sparsity.

  return(list("Kj_comp"=Kj,"Kj_all"=Kj_all))
}
