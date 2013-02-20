
##' Long-wavelength polarisability approximation for metal-coated spheres (shells)
##' 
##' @title polarisability_shell
##' @param radius in nm
##' @param wavelength vector, in nm
##' @param epsilon complex vector,  shell dielectric function
##' @param f filling ratio
##' @param core refractive index of core
##' @param medium incident medium, real
##' @param ... unused
##' @export
##' @family user_level polarisability shell approximate
##' @author baptiste Auguie
polarisability_shell <- function(radius, wavelength, epsilon, f=1, core=1.5, medium=1.5, ...){

  thickness <- radius*(1-f)
  ss <- sqrt(epsilon) / medium
  sc <- core / medium
  
  k <- 2*pi*medium / wavelength
  x <- k*radius
  
  epsa <- sc^2*(1+2*f^3) + 2*ss^2*(1-f^3)
  epsb <- sc^2*(1-f^3) + ss^2*(2+f^3)
  r <- epsa / epsb
  
  epsab <- sc^2*(2+3*f^5) + 3*ss^2*(1-f^5)
  epsbb <- 2*sc^2*(1-f^5) + ss^2*(3+2*f^5)
  rb <- epsab / epsbb

  ## Delta1 <- 2i/3 * x^3* (r*ss^2-1) / 
  ##   (r*ss^2 + 2 - 3/5*x^2 * ((ss^2-1)*((3*r-2)*ss^2-2)+3*f^2*ss^4/epsb*(1-r)*(2*ss^2-sc^2)) /
  ##    (r*ss^2 - 1) - 2i/3 * x^3* (r*ss^2-1))

  ## dip <- 
  ## transform(data.frame(wavelength=wavelength,
  ##                      extinction = -2/x^2 * 3* Re(Delta1),
  ##                      scattering = 2/x^2 * 3* Mod(Delta1)^2),
  ##           absorption = extinction - scattering,
  ##           order = "dipole")   
  
    d10 <-  2i/3 * x^3* (r*ss^2-1) / (r*ss^2+2)
    a1 <- ((ss^2-1)*((3*r-2)*ss^2-2)+3*f^2*ss^4/epsb*(1-r)*(2*ss^2-sc^2)) /
      ((r*ss^2 - 1)*(r*ss^2 +2))
    a2 <- 1 / ((r*s^2 - 1)*(r*s^2+2)) *
      (
       (s^2-1)/2*((43*r^2 - 73*r +32)*s^4 + (25 - 73*r)*s^2 + 32) +
       3*f^4*s^4/epsb*(r-1)*(sc^4 - 24*sc^2*s^2 + 16*s^4) +
       126*(s^2*(s^2-1)*(r-1)) / (r*s^2-1) * (1/6*((4-3*r)*s^2+1)+f^2*s^2/epsb*(sc^2-2*s^2))^2
       )
    a2 <- 1 / ((r*ss^2 - 1)*(r*ss^2+2)) *
      (
       (ss^2-1)/2*((43*r^2 - 73*r +32)*ss^4 + (25 - 73*r)*ss^2 + 32) +
       3*f^4*ss^4/epsb*(r-1)*(sc^4 - 24*sc^2*ss^2 + 16*ss^4) +
       126*(ss^2*(ss^2-1)*(r-1)) / (r*ss^2-1) * (1/6*((4-3*r)*ss^2+1)+f^2*ss^2/epsb*(sc^2-2*ss^2))^2
       )
  Delta14 <- d10 / (1 -  3/5*x^2 * a1 - 3/350*x^4*a2 - d10)
  
  dip <- 
  transform(data.frame(wavelength=material$wavelength,
                       extinction = -2/x^2 * 3* Re(Delta14),
                       scattering = 2/x^2 * 3* Mod(Delta14)^2),
            absorption = extinction - scattering,
            order = "dipole")   
  

  num2 <- 1i/30*x^5*(rb*ss^2-1)
  denom2 <- rb*ss^2 + 3/2 +
    5/14*x^2 * ((ss^2-1)*(1 + (1-rb)*ss^2) + 5*f^2 * ss^6/epsbb*(rb-1)) / 
      (rb*ss^2-1) -1i/30*x^5*(rb*ss^2-1)
  
  Delta2 <- num2 / denom2

  quad <- 
  transform(data.frame(wavelength=wavelength,
                       extinction = -2/x^2 * 5 * Re(Delta2),
                       scattering = 2/x^2 * 5* Mod(Delta2)^2),
            absorption = extinction - scattering,
            order = "quadrupole")
    
  all <- transform(dip,
                   extinction = extinction + quad$extinction, 
                   scattering = scattering + quad$scattering,
                   absorption = absorption + quad$absorption, order="all")

  rbind(dip, quad, all)
  
}

