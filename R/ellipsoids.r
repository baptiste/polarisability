##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param laxis 
##' @param saxis 
##' @param ... unused
##' @return 
##' @export
##' @family user_level polarisability shell approximate
##' @author Baptiste Auguie
chi <- function(laxis=50, saxis=50){
  e <- sqrt(1-saxis^2/laxis^2)
  (1-e^2)/e^2*(-1+1/(2*e)*log((1+e)/(1-e)))
}

##' Long-wavelength polarisability approximation for prolate ellipsoids
##' 
##' @title polarisability_ellipsoid
##' @param wavelength vector, in nm
##' @param epsilon complex vector,  shell dielectric function
##' @param a 
##' @param b 
##' @param medium incident medium, real
##' @param ... unused
##' @export
##' @family user_level polarisability shell approximate
##' @author baptiste Auguie
polarisability_ellipsoid <- function(wavelength, epsilon, a=50, b=30, 
                                     medium = 1.33, ...){
  La <- if(a == b) 1/3 else chi(a, b)
  Lb <- Lc <-  1/2 * (1 - La)
    
  Vcm <- a*b*b
  V <- 4*pi/3*Vcm

  Vcm*(epsilon - medium^2) * 1 / 
    cbind(medium^2 + La*epsilon, medium^2 + Lb*epsilon, medium^2 + Lc*epsilon)
}
