#' @title Simulation of a stock
#' @description \code{stockSim} simulates a stock using variable fishing mortality 
#' and stock recruitment relationship
#' 
#' @param params List of parameters to for the population. Applies a single cohort
#' simulation for the initial population state using \code{\link[fishdynr]{cohortSim}}
#' function
#' @param nyears number of years in the simulation
#' @param Ft time series vector for fishing mortality. If a single value, then the function
#' assumes a constant fishing mortality for the entire simulation (default=0)
#' @param envKt time series vector for environmental effects to maximum recruitment 
#' (e.g. \code{rmax} in \code{\link[fishdynr]{srrBH}}) (default=1).
#' @param envSt time series vector for environmental effects to half maximum recruitment
#' parameter (e.g. \code{beta} in \code{\link[fishdynr]{srrBH}}) (default=1).
#' 
#' @details \code{params} list should contain the following parameters:
#' \itemize{
#'   \item species. Species name
#'   \item growthFun. Name of growth function ("growth_VB" is von Bertalanffy growth function)
#'   \item K. Growth constant (for use in von Bertalanffy growth function))
#'   \item Linf. Infinite length (for use in von Bertalanffy growth function)
#'   \item t0. (hypothetical) time when length equals zero (for use in von Bertalanffy growth function)
#'   \item amax. Maximum age
#'   \item LWa. Length-weight relationship parameter a (weight~a*length^b)
#'   \item LWb. Length-weight relationship parameter b (weight~a*length^b)
#'   \item M. Natural mortality
#'   \item F. Fishing mortality
#'   \item N0. Number of individuals at time 0
#'   \item matFun. Name of maturity function ("pmat_w" is a logistic function that includes width, w, of quantiles)
#'   \item Lmat. Length at maturity (i.e. where probability of being mature is 50%) 
#'   (for use in "pmat_w" function)
#'   \item wmat. Width of length between 25% and 75% maturity. Defines steepness
#'   of transition from immature to mature (for use in "pmat_w" function) 
#'   \item fec. Number of eggs produced per weight [g] of mature female 
#'   \item selectFun. Function to use for gear selection. Determines lengths 
#'   vulnerable to fishing mortality.
#'   \item ... Other parameters for maturity and selectivity functions.
#' }
#' 
#' @return A list
#' 
#' @examples
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
#' tmp <- stockSim(Ft=Ft, params=params, nyears=100, envKt=envKt, envSt=envSt)
#' plot(tmp$Bt, t="l")
#' plot(tmp$Yt, t="l")
#' sum(tmp$Yt/1e6, na.rm=TRUE)
#' 
#' @export
#' 
stockSim <- function(params, nyears=100, Ft=0, envKt=1, envSt=1){
  params$F <- Ft[1]
  res <- cohortSim(params, t_incr=1) # initial stock size
  if(length(Ft)==1) Ft <- rep(Ft, nyears)
  if(length(envKt)==1) envKt <- rep(envKt, nyears)
  if(length(envSt)==1) envSt <- rep(envSt, nyears)
  
  L <- matrix(0, length(res$t), length(res$t)) # Leslie matrix
  subdiag <- which(row(L) == col(L) + 1) # position of subdiagonal
  Ntc <- matrix(0, nrow=nyears, ncol=length(res$t))
  Ntc[1,] <- res$Nt
  Ctc <- matrix(NaN, nrow=nyears, ncol=length(res$t))
  for(i in 2:nyears){
    Neggs <- sum(Ntc[i-1,] * res$pmat * res$Neggst)
    args.incl <- which(names(res) %in% names(formals(get(res$srrFun))))
    Nrecr <- do.call( get(res$srrFun), list(rmax=res$rmax*envKt[i], beta=res$beta*envSt[i], Neggs=Neggs))  
    Ft.i <- Ft[i]*res$pcap
    Zt.i <- res$M + Ft.i
    Mt.i <- Zt.i - Ft.i
    L[subdiag] <- exp(-Zt.i)[-length(Zt.i)] # subdiagonal substitution of survivorship values (gx's)
    Ntc[i,] <- L %*% Ntc[i-1,]
    L[subdiag] <- exp(-Mt.i)[-length(Mt.i)] # subdiagonal substitution of survivorship values (gx's)
    Ntc.noF <- L %*% Ntc[i-1,]
    Ctc[i-1,] <- Ntc.noF - Ntc[i,]  
    Ntc[i,1] <- Nrecr
  }
  
  Btc <- t(apply(Ntc, 1, function(x) x*res$Wt))
  Ytc <- t(apply(Ctc, 1, function(x) x*res$Wt))
  Bt <- rowSums(Btc, na.rm=TRUE)
  Yt <- rowSums(Ytc, na.rm=TRUE)
  Nt <- rowSums(Ntc, na.rm=TRUE)
  Ct <- rowSums(Ctc, na.rm=TRUE)
  
  res2 <- list(
    t=seq(nyears), Ft=Ft, envKt=envKt, envSt=envSt, 
    Btc=Btc, Ntc=Ntc, Ctc=Ctc, Ytc=Ytc,
    Bt=Bt, Nt=Nt, Ct=Ct, Yt=Yt
  ) 
  return(res2)
}
