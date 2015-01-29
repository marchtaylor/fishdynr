#' @title Stock-recruitment relationship (Beverton-Holt type)
#' @description \code{srrBH} describes stock-recruitment relationship 
#' as the number of recruits resulting from a given spawning biomass.
#' 
#' @param rmax maximum recruitment level
#' @param beta parameter describing steepness of relationship. Specifically, 
#' \code{beta} describes the point where the number of spawned eggs results in
#' half the maximum number of recruits, \code{rmax}
#' @param SB spawning biomass of the adult population
#' 
#' @return
#' Number of recruited individuals
#'   
#' @examples
#' SB <- seq(0,1e11,,100)
#' rmax = 2e8
#' beta = 1e10
#' Nrecr <- srrBH(rmax, beta, SB)
#' plot(SB, Nrecr, t="l", ylim=c(0, rmax))
#' abline(h=rmax, lty=2, col=8)
#' lines(x=c(0, beta, beta), y=c(rmax/2, rmax/2, 0), lty=2, col=8)
#' text(x=0, y=rmax*0.95, labels="rmax", col=8, pos=4)
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
srrBH <- function(rmax=1000, beta=500, SB=500){
  nrecr <- rmax*SB/(beta+SB)
  nrecr
}

