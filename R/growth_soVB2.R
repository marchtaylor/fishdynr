#' @title seasonally oscillating von Bertalanffy growth function (using L0)
#' @description \code{growth_soVB2} describes the growth as a function of age (t)
#' using the seasonally oscillating von Bertalanffy growth function (Somers 1988).
#' The function differs slightly from \code{\link[fishdynr]{growth_soVB}} 
#' by using \code{L0} instead of \code{t0} to describe the growth function origin.
#' 
#' @param Linf Infinite length
#' @param K growth constant
#' @param t age 
#' @param L0 (hypothetical) length at time zero
#' @param ts summer point. Time of year (between 0 and 1) when growth oscillation
#' cycle begins (sine wave term becomes positive). Note that this definition 
#' differs from some interpretations of the model (see Somers 1998)
#' @param C oscillation strength. Varies between 0 and 1.
#'   
#' @examples
#' t <- seq(0,5,0.1)
#' L <- growth_soVB2(Linf=100, K=0.5, t=t, L0=10, ts=0.5, C=0.5)
#' plot(t, L, t="l", ylim=c(0,110), xaxs="i", yaxs="i")
#' points(0,10); abline(h=100, lty=2, col=8)
#' 
#' @export
#' 
growth_soVB2 <- function(Linf, K, t, L0, ts, C){
  Linf - (Linf - L0) * exp(-(
    K*(t)
    + (((C*K)/(2*pi))*sin(2*pi*(t-ts)))
  ))
}