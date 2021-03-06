% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_Kj_compress.R
\name{get_Kj_compress}
\alias{get_Kj_compress}
\title{Function to get Kj from compressed X}
\usage{
get_Kj_compress(wavespecsc)
}
\arguments{
\item{wavespecsc}{List of specifications for the wavelet transform of X. Specifically, wavespecsc$Kj, the number of coefs per level, and
wavespecsc$keep, are needed.}
}
\value{
This outputs a list of two values: "Kj_comp" which is # of coefficients per level after the compression, and "Kj_all" which is
# of coefficients per level after the original wavelet transformation.
}
\description{
Given wavespecsc object with Kj = # coeffs per level and keep = vector of 0's and 1's with 1 in coordinate for coefficients to "keep",
figures out J  = # levels and Kj = # coefficients per level for the reduced set of coefficients that are kept after compression.
}
\examples{

}
