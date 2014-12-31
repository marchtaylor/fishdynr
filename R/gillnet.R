#' @title gillnet selection
#' @description \code{gillnet} describes gillnet selection following 
#' Millar and Holst (1997). Possible selectivity distributions include
#' normal (fixed width) and lognormal. [Note: others will be added] 
#' 
#' @param Lt body size
#' @param mesh_size mesh size
#' @param mesh_size1 smallest reference mesh size
#' @param select_dist selectivity type (\code{"normal_fixed", "lognormal"})
#' @param select_p1 selectivity function parameter 1 (see Millar and Holst 1997)
#' @param select_p2 selectivity function parameter 2 (see Millar and Holst 1997)
#'   
#' @examples
#' data(tilapia)
#' tilapia$selectFun="gillnet"
#' mesh_sizes <- c(60, 80, 100, 120)
#' for(i in seq(mesh_sizes)){
#'   tilapia$mesh_size <- mesh_sizes[i]
#'   res <- cohortSim(tilapia, t_incr=0.01)
#'   if(i == 1) plot(pcap ~ Lt, res, t="n")
#'   lines(pcap ~ Lt, res, col=i)
#' }
#' legend("topleft", legend=mesh_sizes, 
#'   col=seq(mesh_sizes), lty=1, 
#'   title="mesh size [mm]", bty="n"
#' )
#' 
#' 
#' @references
#' Millar, R. B., & Holst, R. (1997). 
#' Estimation of gillnet and hook selectivity using log-linear models. 
#' ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export
#' 
gillnet <- function(Lt, mesh_size, mesh_size1, select_dist, select_p1, select_p2){
  if(select_dist == "lognormal"){
    sel <- (1/Lt) * exp(select_p1 + log(mesh_size/mesh_size1) - (select_p2^2/2) - (((log(Lt) - select_p1 - log(mesh_size/mesh_size1))^2) / (2*select_p2^2) ) )
  }
  if(select_dist == "normal_fixed"){
    sel <- exp(-((Lt-mesh_size*select_p1)^2/(2*select_p2^2)))
  } 
  return(sel)
}

