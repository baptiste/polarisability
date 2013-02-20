
##' Long-wavelength polarizability approximation for homogeneous spheres
##' 
##' @title polarizability_sphere
##' @param radius in nm
##' @param wavelength vector, in nm
##' @param epsilon complex vector
##' @param medium incident medium, real
##' @param ... unused
##' @export
##' @family user_level polarizability sphere approximate
##' @author baptiste Auguie
polarizability_sphere <- function(radius, wavelength, epsilon, medium=1.336, ...){

  s <- sqrt(epsilon) / medium
  k <- 2*pi*medium / wavelength
  x <- k*radius
  s2 <- s^2
  D1_0 <- 2i/3 * (s2 - 1)/(s2 + 2) * x^3
  Delta <-  D1_0 / (1 - 3/5*x^2*(s2 - 2)/(s2 + 2) - D1_0 - 3/350*x^4 * (s^4 - 24*s2 + 16) / (s2 + 2))
    
  dip <- 
  transform(data.frame(wavelength=wavelength,
                       extinction = -2/x^2 * 3* Re(Delta),
                       scattering = 2/x^2 * 3* Mod(Delta)^2),
            absorption = extinction - scattering,
            order = "dipole")   

  Delta2 <- (1i*x^5/30*(s2-1)) / (s2+3/2+5/14*x^2-5/2646*x^4*(s^4+30*s2-45)-1i*x^5/30*(s2-1))
  
  quad <- 
  transform(data.frame(wavelength=wavelength,
                       extinction = -2/x^2 * 5 * Re(Delta2),
                       scattering = 2/x^2 * 5* Mod(Delta2)^2),
            absorption = extinction - scattering,
            order = "quadrupole")
  
  Delta3 <- (4i/4725*x^7*(s2-1)) / (s2 +4/3+7*x^2/135*(s2+4) - 7*x^4/10692*(s^4+8*s^2-32) - 4i*x^7/4725*(s2-1))
  oct <- 
  transform(data.frame(wavelength=wavelength,
                       extinction = -2/x^2 * 7 * Re(Delta3),
                       scattering = 2/x^2 * 7 * Mod(Delta3)^2),
            absorption = extinction - scattering,
            order = "octupole")

  
  all <- transform(dip,
                   extinction = extinction + quad$extinction + oct$extinction, 
                   scattering = scattering + quad$scattering + oct$scattering,
                   absorption = absorption + quad$absorption + oct$absorption, order="all")

  rbind(dip, quad, oct, all)
  
}
