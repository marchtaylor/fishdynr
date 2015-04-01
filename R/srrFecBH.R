#' @title Stock-recruitment relationship using fecundity (Beverton-Holt type)
#' @description \code{srrFecBH} describes stock-recruitment relationship 
#' as the number of recruits resulting from a given number of eggs produced by the population.
#' 
#' @param srrFecBH_a long-term mean survival parameter
#' @param srrFecBH_b carrying capacity parameter 
#' @param neggs number of eggs spawned
#' 
#' @return
#' Number of recruited individuals
#'   
#' @examples
#' n <- 100
#' neggs <- seq(0,1e8,,n)
#' srrFecBH_a = 0.2
#' srrFecBH_b = 1e7/srrFecBH_a
#' Nrecr <- srrFecBH(srrFecBH_a, srrFecBH_b, neggs)
#' plot(neggs, Nrecr, t="l", ylim=c(0, srrFecBH_b*srrFecBH_a))
#' abline(0, srrFecBH_a, col=8, lty=2)
#' abline(h=srrFecBH_b*srrFecBH_a, lty=2, col=8)
#' 
#' @references
#' Beverton, R. J., Holt, S. J., 1957. 
#' On the dynamics of exploited fish populations
#'
#' @export
#' 
srrFecBH <- function(srrFecBH_a, srrFecBH_b, neggs){
  nrecr <- (srrFecBH_a*neggs) / (1+(neggs/(srrFecBH_b)))
  nrecr
}
