#' @title Simulation of a cohort
#' @description \code{cohortSim} simulates a single cohort
#' 
#' @param params List of parameters to use in function
#' @param t_incr Value to use a time increment (in years). Defaults to 1.
#' 
#' @description \code{params} list should contain the following parameters:
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
#' res <- cohortSim(tilapia, t_incr=.1)
#' plot(pcap ~ Lt, res, t="l")
#' plot(Lt ~ t, res, t="l")
#' plot(Wt ~ t, res, t="l")
#' 
#' plot(Bt ~ t, res, t="l")
#' lines(SBt ~ t, res, col=2)
#' 
#' plot(Bt ~ Lt, res, t="l")
#' lines(SBt ~ Lt, res, col=2)
#' 
#' plot(Yt ~ t, res, t="l")
#' 
#' plot(Nt ~ t, res, t="l", log="y")
#' lines(Nt.noF ~ t, res, col=2, lty=2)
#' 
#' @export
#' 
cohortSim <- function(params, t_incr=1){
res <- params
if(is.null(res$amax)) res$amax <- with(res, ceiling(log(-(0.95* Linf)/Linf + 1) / -K + t0))
# Ages
t <- seq(0, res$amax, t_incr)
res$t <- t
# Growth  
args.incl <- which(names(res) %in% names(formals(get(res$growthFun))))
Lt <- do.call(get(res$growthFun), res[args.incl])
res$Lt <- Lt
Wt <- with(res, LWa*Lt^LWb)
# Maturity
args.incl <- which(names(res) %in% names(formals(get(res$matFun))))
pmat <- do.call(get(res$matFun), res[args.incl])
Neggst <- Wt * pmat * res$fec/2
# Prob. of capture
args.incl <- which(names(res) %in% names(formals(get(res$selectFun))))
pcap <- do.call(get(res$selectFun), res[args.incl])
#Numbers
N0 <- res$N0; M <- res$M; F <- res$F
Nt.noF <- with(res, N0 * exp(-M*t))
Nt <- 0*t
Nt[1] <- N0
Ct <- 0*t
Nt[1] <- N0
for(i in 2:length(t)){
  Nt[i] <- Nt[i-1] * exp(-(M + F*pcap[i]) * (t[i]-t[i-1]))
  Ct[i-1] <- Nt[i-1] * (1 - exp(-(F*pcap[i]) * (t[i]-t[i-1])))
}
#Yield
Yt <- Ct * Wt
#Population weight
Bt <- Nt * Wt
Bt.noF <- Nt.noF * Wt
#Spawning pop biomass and fecundity
SBt <- Bt*pmat
FECt <- 
  #Optimal length
  Lopt <- Lt[which.max(Bt.noF)]
Lopt.plus_minus_10 <- c(Lopt-Lopt*0.1, Lopt+Lopt*0.1)
Lt.in.Lopt <- Lt >= Lopt.plus_minus_10[1] & Lt <= Lopt.plus_minus_10[2]
Lt.in.mega <- Lt >= Lopt.plus_minus_10[2]
#Stats
SB <- sum(SBt) / t_incr
Y <- sum(Yt, na.rm=TRUE)
Ftot <- sum(Yt/Bt, na.rm=TRUE) / res$amax  
fracC.mat <- sum(Ct * pmat, na.rm=TRUE) / sum(Ct, na.rm=TRUE)
fracC.Lopt <- sum(Ct * Lt.in.Lopt, na.rm=TRUE) / sum(Ct, na.rm=TRUE)
fracN.mega <- sum(Nt * Lt.in.mega, na.rm=TRUE) / sum(Nt, na.rm=TRUE)
fracC.mega <- sum(Ct * Lt.in.mega, na.rm=TRUE) / sum(Ct, na.rm=TRUE)
# Out
res2 <- list(
  Wt=Wt,
  Nt=Nt, Nt.noF=Nt.noF, Bt=Bt,
  pmat=pmat, Neggst=Neggst, pcap=pcap,
  Ct=Ct, SBt=SBt, Yt=Yt,
  Lopt=Lopt,
  SB=SB, Y=Y,
  fracC.mat=fracC.mat, fracC.Lopt=fracC.Lopt,
  fracN.mega=fracN.mega, fracC.mega=fracC.mega
)
return(c(res, res2))  
}