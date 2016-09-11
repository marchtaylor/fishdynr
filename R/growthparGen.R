#' Generate growth parameters for soVBGF
#' 
#' @description Generate growth parameters for a population for the seasonally oscillating 
#' von Bertalanffy growth function with a given given mean and coefficient 
#' of variation. Generated parameters can be passed to \code{\link[fishdynr]{growth_soVB}}.
#'
#' @param n number of individuals to general 
#' @param K.mu mean growth constant K of population
#' @param K.cv coefficient of variation for growth constant K
#' @param Linf.mu mean infinite length of population
#' @param Linf.cv coefficient of variation for infinite length Linf
#' @param ts summer point (0-1)
#' @param C oscilation strength (0-1)
#' @param t0 age at lengthh zero
#'
#' @return data.frame with seasonally oscillating von Bertalanffy growth function 
#' parameters for each individual (rows).
#' 
#' @references
#' Vakily, J.M., 1992. Determination and comparison of bivalve growth, 
#' with emphasis on Thailand and other tropical areas. WorldFish.
#' 
#' Munro, J.L., Pauly, D., 1983. A simple method for comparing the growth 
#' of fishes and invertebrates. Fishbyte 1, 5-6.
#' 
#' Pauly, D., Munro, J., 1984. Once more on the comparison of growth 
#' in fish and invertebrates. Fishbyte (Philippines).
#' 
#' @export
#'
#' @examples
#'
#' set.seed(1111)
#' inds <- growthparGen(n=500)
#' 
#' # check of growth performance index (phi')
#' plot(K ~ Linf, inds, log="xy")
#' coef(lm(log10(K) ~ log10(Linf), inds))[1] # phiprime
#' mean(inds$phiprime) # comparison
#' 
#' # generate growth curves
#' z <- inds$Linf
#' pal <- colorRampPalette(c("#0000FF80", "#00FFFF80", "#FFFF0080", "#FF000080"), alpha=TRUE )
#' col <- pal(100)
#' zlim <- range(z)
#' breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
#' CUT <- cut(z, breaks=breaks, include.lowest = TRUE)
#' colorlevels <- col[match(CUT, levels(CUT))]
#' t <- seq(-0.1,6,len=100)
#' for(i in seq(nrow(inds))){
#'   tmp <- as.list(inds[i,c("Linf", "K", "t0", "ts", "C")])
#'   tmp$t <- t
#'   tmp$Lt <- do.call(growth_soVB, tmp)
#'   if(i == 1){
#'     plot(Lt ~ t, data=tmp, t="n", ylab="Lt", ylim=c(0,max(zlim)))
#'     usr <- par()$usr
#'     rect(usr[1], usr[3], usr[2], usr[4], col="grey85")
#'     grid(col="white")
#'   }
#'   lines(Lt ~ t, data=tmp, col=colorlevels[i])
#' }
#' 
#' 
growthparGen <- function(
n = 100,
K.mu = 0.5, K.cv = 0.1,
Linf.mu = 80, Linf.cv = 0.1,
ts = 0.25, C = 0.85,
t0 = -0.1
){
  
inds <- data.frame(
  Linf = Linf.mu * rlnorm(n, 0, Linf.cv)
)

# mean phiprime
phiprime.mu = log10(K.mu) + 2*log10(Linf.mu)


inds$K <- 10^(phiprime.mu - 2*log10(inds$Linf)) * rlnorm(n, 0, K.cv)
inds$phiprime <- log10(inds$K) + 2*log10(inds$Linf)
inds$t0 <- t0
inds$ts <- ts
inds$C <- C
return(inds)

}
