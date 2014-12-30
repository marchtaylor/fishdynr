#' @title Optimization of fishing mortality
#' @description \code{optim.stockSim} Sets up an optimization for fishing
#'   mortality policy
#'   
#' @param Ft time series vector for fishing mortality. If a single value, then
#'   the function assumes a constant fishing mortality for the entire simulation
#'   (default=0)
#' @param params List of parameters to for the population. Applies a single
#'   cohort simulation for the initial population state using
#'   \code{\link[fishdynr]{cohortSim}} function. See
#'   \code{\link[fishdynr]{stockSim}} for details.
#' @param nyears number of years in the simulation
#' @param envKt time series vector for environmental effects to maximum
#'   recruitment (e.g. \code{rmax} in \code{\link[fishdynr]{srrBH}})
#'   (default=1).
#' @param envSt time series vector for environmental effects to half maximum
#'   recruitment parameter (e.g. \code{beta} in \code{\link[fishdynr]{srrBH}})
#'   (default=1).
#'   
#' @examples
#' \donttest{
#' data(tilapia)
#' params <- tilapia
#' params$knife_edge_size <- 20
#' params$rmax <- 1e6
#' params$beta <- 1e9
#' params$fec <- 100
#' params$N0 <- 1e5
#' nyears <- 100
#' Ft <- rep(0.5, nyears); Ft[20:40] <- 1
#' envKt <- rep(1, nyears); envKt[50:100] <- 0.5
#' envSt <- runif(nyears, min=0.8, max=1.2)
#' 
#' # Optimization of Ft
#' out <- optim(
#' par = Ft,  # initial guess
#' fn = optim.stockSim,
#' params = params,
#' nyears = nyears,
#' envKt = envKt,
#' envSt = envSt,
#' method = "L-BFGS-B",
#' lower = 0,
#' upper = 3,
#' control = list(fnscale=-1, trace=4)
#' )
#' 
#' plot(out$par, t="l")
#' 
#' tmp <- stockSim(params, nyears=100, Ft=out$par, envKt=envKt, envSt=envSt)
#' plot(Yt ~ t, tmp, t="l")
#' sum(tmp$Yt/1e6, na.rm=TRUE)
#' 
#' plot(Yt ~ Bt, tmp)
#' fit <- lm(Yt ~ Bt, tmp)
#' abline(fit)
#' summary(fit)
#' xint <- -coef(fit)[1]/coef(fit)[2]
#' plot(Bt ~ t, tmp, t="l", ylim=c(0, max(tmp$Bt)))
#' lines(Yt ~ t, tmp, col=2)
#' }
#' 
#' @export
#' 
optim.stockSim <- function(Ft=0, params, nyears=100, envKt=1, envSt=1){
  res <- stockSim(params=params, nyears=nyears, Ft=Ft, envKt=envKt, envSt=envSt)
  sum(res$Yt, na.rm=TRUE)
}
