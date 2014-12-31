#' @title knife-edge selection
#' @description \code{knife_edge} describes knife-edge selection where 
#' probability of capture 100\code{\%} after a minimum defined size. Provides
#' a rough estimate of trawl-type selectivity.
#'  
#' @param Lt body size
#' @param knife_edge_size knife edge size. Minimum size at capture.
#'    
#' @examples
#' data(tilapia)
#' tilapia$selectFun="knife_edge"
#' knife_edge_sizes <- c(20, 25, 30, 35)
#' for(i in seq(knife_edge_sizes)){
#'   tilapia$knife_edge_size <- knife_edge_sizes[i]
#'   res <- cohortSim(tilapia, t_incr=0.01)
#'   if(i == 1) plot(pcap ~ Lt, res, t="n")
#'   lines(pcap ~ Lt, res, col=i)
#' }
#' legend("topleft", legend=knife_edge_sizes, 
#'   col=seq(knife_edge_sizes), lty=1, 
#'   title="min. size", bty="n"
#' )
#' 
#' @export
#' 
knife_edge <- function (Lt, knife_edge_size) 
{
  sel <- rep(0, length(Lt))
  sel[Lt >= knife_edge_size] <- 1
  return(sel)
}