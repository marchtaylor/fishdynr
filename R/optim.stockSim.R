#' @title Optimization of fishing mortality
#' @description Sets up an exploration of optimum fishing policy using the function
#' \code{\link[fishdynr]{stockSim}}
#' (i.e. optimal time series of fishing mortalities, \code{Ft}, in order to maximize
#' time series yield, \code{Yt}). The routine is described in further detail in 
#' Walters and Martell (2004).
#'   
#' @param Ft time series vector for fishing mortality. If a single value, then
#'   the function assumes a constant fishing mortality for the entire simulation
#'   (default=0)
#' @param params List of parameters to for the population. Applies a single
#'   cohort simulation for the initial population state using
#'   \code{\link[fishdynr]{cohortSim}} function. See
#'   \code{\link[fishdynr]{stockSim}} for details.
#' @param nyears number of years in the simulation
#' @param env_at time series vector for environmental effects to maximum
#'   recruitment (e.g. \code{srrFecBH_a} in \code{\link[fishdynr]{srrFecBH}})
#'   (default=1).
#' @param env_bt time series vector for environmental effects to maximum
#'   recruitment (e.g. \code{srrFecBH_b} in \code{\link[fishdynr]{srrFecBH}})
#'   (default=1).
#' @param opt what should be optimized ("sum", "sum.log", "sum.disc")
#' @param disc.rate discount rate (default = 0.02)
#'   
#' @references
#' Walters, C. J., Martell, S. J., 2004. 
#' Fisheries ecology and management. Princeton University Press.
#'   
#' @examples
#' \donttest{
#' set.seed(1)
#' data(tilapia)
#' params <- tilapia
#' params$knife_edge_size <- 20
#' params$N0 <- 1e9
#' nyears <- 50
#' Ft <- rep(0.5, nyears)
#' env_at <- runif(nyears, min=0.5, max=1.5)
#' env_bt <- rep(1, nyears); env_bt[20:35] <- 0.5
#' 
#' 
#' # Optimization of Ft (will take some time to reach cost function mimimum)
#' out <- optim(
#'   par = Ft,  # initial guess
#'   fn = optim.stockSim,
#'   params = params,
#'   nyears = nyears,
#'   env_at = env_at,
#'   env_bt = env_bt,
#'   opt = "sum.log", # optimize the sum of log yield
#'   method = "L-BFGS-B",
#'   lower = 0,
#'   upper = 2,
#'   control = list(fnscale=-1, trace=4)
#' )
#' 
#' # optimum Ft series
#' plot(out$par, t="l")
#' 
#' # optimum Yt series
#' tmp <- stockSim(params, nyears=nyears, Ft=out$par, env_at=env_at, env_bt=env_bt)
#' plot(Yt ~ t, tmp, t="l")
#' sum(tmp$Yt/1e6, na.rm=TRUE)
#' 
#' # optimum yield versus stock biomass
#' plot(Yt ~ Bt, tmp)
#' fit <- lm(Yt ~ Bt, tmp)
#' abline(fit)
#' summary(fit)
#' plot(Bt ~ t, tmp, t="l", ylim=c(0, max(tmp$Bt)))
#' lines(Yt ~ t, tmp, col=2)
#' }
#' 
#' @export
#' 
optim.stockSim <- function(Ft=0, params, nyears=100, env_at=1, env_bt=1, opt="sum", disc.rate=0.02){
  res <- stockSim(params=params, nyears=nyears, Ft=Ft, env_at=env_at, env_bt=env_bt)
  incl <- which(res$Yt>0)
  if(opt == "sum") tmp <- sum((res$Yt[incl]), na.rm=TRUE)
  if(opt == "sum.log") tmp <- sum(log(res$Yt[incl]), na.rm=TRUE)
  if(opt == "sum.disc") tmp <- sum((res$Yt / ((1+disc.rate)^(nyears:1)))[incl], na.rm=TRUE)
  return(tmp)
}
