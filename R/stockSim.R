#' @title Simulation of a stock
#' @description \code{stockSim} simulates a stock using variable fishing mortality 
#' and stock recruitment relationship
#' 
#' @param params List of parameters for the population. Applies a single cohort
#' simulation for the initial population state using \code{\link[fishdynr]{cohortSim}}
#' function
#' @param nyears number of years in the simulation
#' @param Ft time series vector for fishing mortality. If a single value, then the function
#' assumes a constant fishing mortality for the entire simulation (default=0)
#' @param env_at time series vector for environmental effects to maximum recruitment 
#' (e.g. \code{srrFecBH_a} in \code{\link[fishdynr]{srrFecBH}}) (default=1).
#' @param env_bt time series vector for environmental effects to half maximum recruitment
#' parameter (e.g. \code{srrFecBH_b} in \code{\link[fishdynr]{srrFecBH}}) (default=1).
#' 
#' @details \code{params} list should contain the following parameters:
#' \itemize{
#'   \item \code{species} Species name
#'   \item \code{growthFun} Name of growth function (e.g. "growth_VB" is the von 
#'   Bertalanffy growth function)
#'   \item \code{amax} Maximum age
#'   \item \code{LWa} Length-weight relationship parameter a (weight~a*length^b)
#'   \item \code{LWb} Length-weight relationship parameter b (weight~a*length^b)
#'   \item \code{M} Natural mortality
#'   \item \code{F} Fishing mortality
#'   \item \code{N0} Number of individuals at time 0
#'   \item \code{matFun} Name of maturity function (e.g. "pmat_w" is a logistic 
#'   function that includes width, w, of quantiles)
#'   \item \code{selectFun} Function to use for gear selection. Determines lengths 
#'   vulnerable to fishing mortality (e.g. "gillnet" and "knife_edge" functions).
#'   \item \code{srrFun} Stock-recruitment relationship function (e.g. "srrFecBH").
#'   \item \code{fec} Number of eggs produced per weight [g] of mature female (For use in
#'   srrFecFun).
#'   \item \code{...} Other parameters for growth, maturity, and selectivity functions.
#' }
#' For fitting an optimal time series of fishing mortalities, \code{Ft}, see
#' \code{\link[fishdynr]{optim.stockSim}} (Walters and Martell, 2004).
#' 
#' 
#' @return A list
#' \itemize{
#'   \item \code{Btc} matrix. Stock biomass by time (rows) and cohort (columns).
#'   \item \code{Ytc} matrix. Fishery yield by time (rows) and cohort (columns).
#'   \item \code{Bt} vector. Stock biomass by time.
#'   \item \code{Yt} vector. Fishery yield by time.
#'   \item \code{Nt} vector. Stock size (in numbers) by time.
#'   \item \code{Ct} vector. Fishery catch (in numbers) by time.
#' }
#' 
#' @references
#' Walters, C. J., Martell, S. J., 2004. 
#' Fisheries ecology and management. Princeton University Press.
#' 
#' @examples
#' data(tilapia)
#' params <- tilapia
#' params$knife_edge_size <- 20
#' params$N0 <- 1e9
#' nyears <- 50
#' Ft <- rep(0.5, nyears)
#' env_at <- runif(nyears, min=0.5, max=1.5)
#' env_bt <- rep(1, nyears); env_bt[20:35] <- 0.5
#' tmp <- stockSim(Ft=Ft, params=params, nyears=nyears, env_at=env_at, env_bt=env_bt)
#' plot(tmp$Bt, t="l")
#' plot(tmp$Yt, t="l")
#' sum(tmp$Yt/1e6, na.rm=TRUE)
#' 
#' @export
#' 
stockSim <- function(params, nyears=100, Ft=0, env_at=1, env_bt=1){
  params$F <- Ft[1]
  res <- cohortSim(params, t_incr=1) # initial stock size
  if(length(Ft)==1) Ft <- rep(Ft, nyears)
  if(length(env_at)==1) env_at <- rep(env_at, nyears)
  if(length(env_bt)==1) env_bt <- rep(env_bt, nyears)
  
  L <- matrix(0, length(res$t), length(res$t)) # Leslie matrix
  subdiag <- which(row(L) == col(L) + 1) # position of subdiagonal
  Ntc <- matrix(0, nrow=nyears, ncol=length(res$t))
  Ntc[1,] <- res$Nt
  Btc <- matrix(0, nrow=nyears, ncol=length(res$t))
  Btc[1,] <- res$Bt
  SBtc <- matrix(0, nrow=nyears, ncol=length(res$t))
  SBtc[1,] <- res$SBt
  Fectc <- matrix(0, nrow=nyears, ncol=length(res$t))
  Fectc[1,] <- res$Neggst * res$Nt * res$pmat * res$nspawn * res$p_female
  Ctc <- matrix(NaN, nrow=nyears, ncol=length(res$t))
  for(i in 2:nyears){
    Nrecr <- do.call(get(res$srrFun), args=list(srrFecBH_a=params$srrFecBH_a*env_at[i], srrFecBH_b=params$srrFecBH_b*env_bt[i], neggs=sum(Fectc[i-1,])))  
    Ft.i <- Ft[i]*res$pcap
    Zt.i <- res$M + Ft.i
    Mt.i <- Zt.i - Ft.i
    L[subdiag] <- exp(-Zt.i)[-length(Zt.i)] # subdiagonal substitution of survivorship values (gx's)
    Ntc[i,] <- L %*% Ntc[i-1,]
    L[subdiag] <- exp(-Mt.i)[-length(Mt.i)] # subdiagonal substitution of survivorship values (gx's)
    Ntc.noF <- L %*% Ntc[i-1,]
    Ctc[i-1,] <- Ntc.noF - Ntc[i,]  
    Ntc[i,1] <- Nrecr
    Btc[i,] <- Ntc[i,] * res$Wt
    SBtc[i,] <- Btc[i,] * res$pmat
    Fectc[i,] <- res$Neggst * Ntc[i,] * res$pmat * res$nspawn * res$p_female
  }
  
  Ytc <- t(apply(Ctc, 1, function(x) x*res$Wt))
  Bt <- rowSums(Btc, na.rm=TRUE)
  SBt <- rowSums(SBtc, na.rm=TRUE)
  Yt <- rowSums(Ytc, na.rm=TRUE)
  Nt <- rowSums(Ntc, na.rm=TRUE)
  Fect <- rowSums(Fectc, na.rm=TRUE)
  Ct <- rowSums(Ctc, na.rm=TRUE)
  
  res2 <- list(
    t=seq(nyears), Ft=Ft, env_at=env_at, env_bt=env_bt, 
    Btc=Btc, SBtc=SBtc, Ntc=Ntc, Fectc=Fectc, Ctc=Ctc, Ytc=Ytc,
    Bt=Bt, SBt=SBt, Nt=Nt, Fect=Fect, Ct=Ct, Yt=Yt
  ) 
  return(res2)
}
