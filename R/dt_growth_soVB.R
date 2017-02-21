#' @title Discrete seasonally oscillating von Bertalanffy growth function
#' @description \code{dt_growth_soVB} calculates length at time 2, given 
#' length at time 1 according to the seasonally oscillationg 
#' von Bertalanffy growth function
#' 
#' @param Linf Infinite length
#' @param K growth constant
#' @param ts summer point. Time of year (between 0 and 1) when growth oscillation
#' cycle begins (sine wave term becomes positive). Note that this definition 
#' differs from some interpretations of the model (see Somers 1998)
#' @param C oscillation strength. Varies between 0 and 1
#' @param L1 Length at t1
#' @param t1 time 1
#' @param t2 time 2
#' 
#'   
#' @examples
#' # Growth parameters
#' Linf = 100
#' K = 0.5
#' ts = 0.5
#' C = 0.75
#' 
#' # create dataframe
#' df <- data.frame(t = seq(0,5,0.2), L=NaN)
#' df$L[1] <- 10 # intial size
#' 
#' # run iterative model
#' for(i in 2:nrow(df)){
#'   df$L[i] <- dt_growth_soVB(
#'     Linf, K, ts, C, 
#'     L1=df$L[i-1],
#'     t1=df$t[i-1],
#'     t2=df$t[i]
#'   )
#' }
#' 
#' # plot result
#' plot(L ~ t, df)
#' arrows(
#'   x0 = df$t[-nrow(df)], y0 = df$L[-nrow(df)],
#'   x1 = df$t[-1], y1 = df$L[-1],
#'   col = 8, length = 0.15
#' )
#' 
#' @references
#' Somers, I. F. (1988). On a seasonally oscillating growth function. 
#' Fishbyte, 6(1), 8-11.
#'
#' @export
#' 
dt_growth_soVB <- function(Linf, K, ts, C, L1, t1, t2){
  dt <- (Linf - L1) *
  {1 - exp(-(
    K*(t2-t1)
    - (((C*K)/(2*pi))*sin(2*pi*(t1-ts)))
    + (((C*K)/(2*pi))*sin(2*pi*(t2-ts)))
  ))}
  L2 <- L1 + dt
  L2
}