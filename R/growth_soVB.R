#' @title seasonally oscillating von Bertalanffy growth function
#' @description \code{growth_soVB} describes the growth as a function of age (t)
#' using the seasonally oscillating von Bertalanffy growth function (Somers 1988).
#' 
#' @param Linf Infinite length
#' @param K growth constant
#' @param t age 
#' @param t0 (hypothetical) age at length zero
#' @param ts summer point. Time of year (between 0 and 1) when growth oscillation
#' cycle begins (sine wave term becomes positive). Note that this definition 
#' differs from some interpretations of the model (see Somers 1998)
#' @param C oscillation strength. Varies between 0 and 1.
#' 
#'   
#' @examples
#' t <- seq(0,5,0.1)
#' L <- growth_soVB(Linf=100, K=0.5, t=t, t0=-0.2, ts=0.5, C=0.75)
#' plot(t, L, t="l")
#' 
#' @references
#' Somers, I. F. (1988). On a seasonally oscillating growth function. 
#' Fishbyte, 6(1), 8-11.
#'
#' @export
#' 
growth_soVB <- function(Linf, K, t, t0, ts, C){
  Linf * (1-exp(-(
    K*(t-t0)
    + (((C*K)/(2*pi))*sin(2*pi*(t-ts)))
    - (((C*K)/(2*pi))*sin(2*pi*(t0-ts)))
  )))
}