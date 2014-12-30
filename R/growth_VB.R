#' @title von Bertalanffy growth function
#' @description \code{growth_VB} describes the growth as a function of age (t)
#' using the von Bertalanffy growth function
#' 
#' @param Linf Infinite length
#' @param K growth constant
#' @param t age 
#' @param t0 (hypothetical) age at length zero
#'   
#' @examples
#' t <- seq(0,5,0.1)
#' L <- growth_VB(Linf=100, K=0.5, t=t, t0=-0.2)
#' plot(t, L, t="l")
#' 
#' @export
#' 
growth_VB <- function(Linf, K, t, t0){
  Linf * (1-exp(-K*(t-t0)))
}