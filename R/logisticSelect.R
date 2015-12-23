#' @title logistic selection
#' @description \code{logisticSelect} describes knife-edge selection where 
#' probability of capture 100\code{\%} after a minimum defined size. Provides
#' a rough estimate of trawl-type selectivity.
#'  
#' @param Lt body size
#' @param L50 length at first capture (i.e. where prob. of capture equals 50\%)
#' @param wqs width of quartiles between 25\% and 75\% probability of capture
#'    
#' @examples
#' data(tilapia)
#' tilapia$selectFun="logisticSelect"
#' tilapia$L50 <- 20
#' WQS <- c(0, 3, 6)
#' COL <- rainbow(3)
#' for(i in seq(WQS)){
#'   tilapia$wqs <- WQS[i]
#'   res <- cohortSim(tilapia, t_incr=0.01)
#'   if(i == 1) plot(pcap ~ Lt, res, t="n")
#'   lines(pcap ~ Lt, res, col=COL[i])
#' }
#' legend("bottomright", 
#'        legend=paste(WQS, c(" ('knife-edge')", "", ""), sep=""), 
#'        col=COL, lty=1, 
#'        title="Width between 0.25 and 0.75 quantiles", bty="n"
#' )
#' lines(c(0,tilapia$L50,tilapia$L50), c(0.5,0.5,0), col=1, lty=2)
#' text(20,.5,bquote("L"[50]), pos=4, col=1, font=1)
#' 
#' 
#' @export
#' 
logisticSelect <- function(Lt, L50, wqs){
  1 / (1 + exp(-(Lt - L50) / (wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))) ) 
}