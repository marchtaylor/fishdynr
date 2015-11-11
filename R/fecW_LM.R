#' @title Fecundity and function of weight
#' @description \code{fecW_LM} is a linear model describing fecundity (number of eggs) as a function
#' of individual weight
#' 
#' @param fec_a intercept
#' @param fec_b slope
#' @param Wt individual weight (from \code{cohortSim} function)
#' 
#' @return
#' Number of eggs per individual
#'   
#' @examples
#' data(tilapia)
#' res <- cohortSim(tilapia, t_incr=0.1)
#' plot(Neggst ~ Wt, res) # Number of eggs as a function of ind. weight
#' plot(res$t, res$Neggst*res$pmat*res$Nt) # total population fecundity by age
#'
#' @export
#' 
fecW_LM <- function(fec_a=0, fec_b=1, Wt=500){
  Neggs <- fec_a + fec_b*Wt
  repl <- which(Neggs < 0)
  if(length(repl) > 0) Neggs[repl] <- 0
  Neggs
}

