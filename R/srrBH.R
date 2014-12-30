#' @title Stock-recruitment relationship (Beverton-Holt type)
#' @description \code{srrBH} describes stock-recruitment relationship
#' 
#' @param rmax maximum recruitment level
#' @param beta parameter describing steepness of relationship. Specifically, 
#' \code{beta} describes the point where the number of spawned eggs results in
#' half the maximum number of recruits, \code{rmax}
#' @param Neggs number of eggs spawned by the adult population
#'   
#' @examples
#' Neggs <- seq(0,5e3,,100)
#' rmax <- 2000
#' beta <- 500
#' Nrecr <- srrBH(rmax, beta, Neggs)
#' plot(Neggs, Nrecr, t="l", ylim=c(0, rmax))
#' abline(h=rmax, lty=2, col=8)
#' lines(x=c(0, beta, beta), y=c(rmax/2, rmax/2, 0), lty=2, col=8)
#' text(x=100, y=rmax*0.95, labels="rmax", col=8)
#' text(x=beta, y=rmax/2, labels="rmax/2", pos=4, col=8)
#' text(x=beta, y=0, labels="beta", pos=4, col=8)
#' 
#' 
#' @references
#' Beverton, R. J., Holt, S. J., 1957. 
#' On the dynamics of exploited fish populations
#'
#' @export
#' 
srrBH <- function(rmax=1000, beta=500, Neggs=500){
  nrecr <- rmax*Neggs/(beta+Neggs) # Michaelis-Menten type
  nrecr
}

