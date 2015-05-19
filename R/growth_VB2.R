#' @title von Bertalanffy growth function (using L0)
#' @description \code{growth_VB2} describes the growth as a function of age (t)
#' using the von Bertalanffy growth function. The function differs slightly from
#' \code{\link[fishdynr]{growth_VB}} by using \code{L0} instead of 
#' \code{t0} to describe the growth function origin
#' 
#' @param Linf Infinite length
#' @param K growth constant
#' @param t age 
#' @param L0 (hypothetical) length at time zero
#'   
#' @examples
#' t <- seq(0,5,0.1)
#' L <- growth_VB2(Linf=100, K=0.5, t=t, L0=10)
#' plot(t, L, t="l", ylim=c(0,110), xaxs="i", yaxs="i")
#' points(0,10); abline(h=100, lty=2, col=8)
#' 
#' @export
#' 
growth_VB2 <- function(Linf, K, t, L0){
  Linf - (Linf - L0)*exp(-K*t)
}