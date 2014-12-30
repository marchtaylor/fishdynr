#' @title Probability of maturity logistic function using quartile width
#' @description \code{pmat_w} describes the probability of maturity as a function 
#' of the size when 50\code{\%}  of the individuals in a population are mature and the width
#' (in size) between 25\code{\%} and 75\code{\%} probability of maturity quartiles 
#' (after Heino et al. 2002).
#' 
#' @param Lt size for probability of maturity calculation
#' @param Lmat size at 50\code{\%} probability of maturity (i.e. "massive maturity")
#' @param wmat width (in size) between 25\code{\%} and 75\code{\%}
#' probability of maturity quartiles
#'   
#' @examples
#' L <- seq(1,20,0.1)
#' pmat1 <- pmat_w(L, Lmat=10, wmat=5)
#' pmat2 <- pmat_w(L, Lmat=10, wmat=2)
#' plot(L, pmat1, t="l", ylab="prob. of maturity")
#' lines(L, pmat2, lty=2)
#' legend("bottomright", legend=c("Lmat=10; wmat=5", "Lmat=10; wmat=2"), 
#'  col=1, lty=1:2, bty="n")
#' 
#' @references
#' Heino M, Dieckmann U, Godo OR (2002) Measuring probabilistic reaction norms 
#' for age and size at maturation. Evolution 56: 669-678.
#'
#' @export
#' 
pmat_w <- function(Lt, Lmat, wmat){
  1 / (1 + exp(-(Lt - Lmat) / (wmat / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))) ) 
}
