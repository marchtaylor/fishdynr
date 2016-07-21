#' @title Generate length frequency data from a synthetic fish population
#'
#' @param tincr time increment
#' @param K.mu mean K (growth parameter from von Bertalanffy growth function) 
#' @param K.cv coefficient of variation on K
#' @param Linf.mu mean Linf (infinite length parameter from von Bertalanffy growth function)
#' @param Linf.cv coefficient of variation on Linf
#' @param ts summer point (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param C strength of seasonal oscillation (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param LWa length-weight relationship constant 'a' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param LWb length-weight relationship constant 'b' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param Lmat length at maturity (where 50\% of individuals are mature)
#' @param wmat width between 25\% and 75\% quantiles for Lmat
#' @param rmax parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param beta parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param repro_wt weight of reproduction (vector of monthly reproduction weight) 
#' @param M natural mortality
#' @param harvest_rate Fishing mortality (i.e. 'F')
#' @param L50 minimum length of capture (in cm). Where selectivity equals 0.5. Assumes logistic ogive typical of trawl net selectivity.
#' @param wqs width of selectivity ogive (in cm)
#' @param bin.size resulting bin size for length frequencies (in cm)
#' @param timemin time at start of simulation (in years). Typically set to zero.
#' @param timemax time at end of simulation (in years).
#' @param timemin.date date corresponding to timemin (of "Date" class)
#' @param tincr time increment for simulation (default = 1/12; i.e. 1 month)
#' @param N0 starting number of individuals
#' @param fished_t times when stock is fished
#' @param lfqFrac fraction of fished stock that are sampled for length frequency data (default = 1).
#' 
#' @description See (see \code{\link[fishdynr]{dt_growth_soVB}})
#' 
#' @return a list containing growth parameters and length frequency object
#' 
#' @export
#'
#' @examples
#' \donttest{
#' res <- lfqGen()
#' names(res)
#' 
#' op <- par(mfcol=c(2,1), mar=c(4,4,1,1))
#' plot(N ~ dates, data=res$pop, t="l")
#' plot(B ~ dates, data=res$pop, t="l", ylab="B, SSB")
#' lines(SSB ~ dates, data=res$pop, t="l", lty=2)
#' par(op)
#' 
#' pal <- colorRampPalette(c("grey30",5,7,2), bias=2)
#' with(res$lfqbin, image(x=dates, y=midLengths, z=t(catch), col=pal(100)))
#' }
#' 
#' 
lfqGen <- function(
tincr = 1/12,
K.mu = 0.5, K.cv = 0.1,
Linf.mu = 80, Linf.cv = 0.1, 
ts = 0, C = 0.85,
LWa = 0.01, LWb = 3,
Lmat = 0.5*Linf.mu, wmat = Lmat*0.2,
rmax = 10000, beta = 1,
repro_wt = c(0,0,1,0,0,0,0,0,0,0,0,0),
M = 0.7, harvest_rate = M, 
L50 = 0.25*Linf.mu, wqs = L50*0.2,
bin.size = 1,
timemin = 0, timemax = 5, timemin.date = as.Date("1980-01-01"),
N0 = 10000,
fished_t = seq(0,5,tincr),
lfqFrac = 1
){

# times
timeseq = seq(from=timemin, to=timemax, by=tincr)
if(!zapsmall(1/tincr) == length(repro_wt)) stop("length of repro_wt must equal the number of tincr in one year")
repro_wt <- repro_wt/sum(repro_wt)
repro_t <- rep(repro_wt, length=length(timeseq)) 
# repro_t <- seq(timemin+repro_toy, timemax+repro_toy, by=1)

# make empty lfq object
lfq <- vector(mode="list", length(timeseq))
names(lfq) <- timeseq

# Estimate tmaxrecr
tmaxrecr <- which.max(repro_wt)*tincr

# Winf.mu and phi
Winf.mu <- LWa*Linf.mu^LWb
phi <- log10(K.mu) + 0.67*log10(Winf.mu)
phiprime = log10(K.mu) + 2*log10(Linf.mu)



# required functions ------------------------------------------------------
date2yeardec <- function(date){as.POSIXlt(date)$year+1900 + (as.POSIXlt(date)$yday)/365}
yeardec2date <- function(yeardec){as.Date(strptime(paste(yeardec%/%1, ceiling(yeardec%%1*365+1), sep="-"), format="%Y-%j"))}

make.inds <- function(
	id=NaN, A = 1, L = 0, W=NaN, mat=0,
	K = K.mu, Winf=NaN, Linf=NaN, phiprime=NaN,
	F=NaN, Z=NaN, 
	Fd=0, alive=1
){
  inds <- data.frame(
    id = id,
    A = A,
    L = L,
    W = W,
    mat = mat,
    K = K,
    Linf = Linf,
    Winf = Winf,
    phiprime = phiprime,
    F = F,
    Z = Z,
    Fd = Fd,
    alive = alive
  )
  lastID <<- max(inds$id)
  return(inds)
}


express.inds <- function(inds){
  inds$Linf <- Linf.mu * rlnorm(nrow(inds), 0, Linf.cv)
  inds$Winf <- LWa*inds$Linf^LWb
  inds$K <- 10^((phi - 0.67*log10(inds$Winf)) * rlnorm(nrow(inds), 0, K.cv))
  #inds$L <- inds$L * rlnorm(nrow(inds), meanlog = 0, sdlog = L0.cv)
  inds$W <- LWa*inds$L^LWb
  inds$phiprime <- log10(inds$K) + 2*log10(inds$Linf)
	return(inds)
}



grow.inds <- function(inds){
	#grow
  L2 <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, C = C, L1 = inds$L, t1 = t-tincr, t2 = t)
  #update length and weight
	inds$L <- L2
	inds$W <- LWa*inds$L^LWb
	#age inds
	inds$A <- inds$A + tincr
	return(inds)
}

mature.inds <- function(inds){
	p <- pmat_w(inds$L, Lmat, wmat) # probability of being mature at length L
	p0t <- (1-p)^(tincr) # prob. of no maturity in tincr
	p1t <- 1-p0t  # probability of maturity in tincr
	inds$mat <- ifelse(runif(nrow(inds)) < p1t | inds$mat, 1, 0)
	return(inds)
}

death.inds <- function(inds){
  pSel <- logisticSelect(inds$L, L50, wqs)
  inds$F <- pSel * Fmax
  inds$Z <- M + inds$F
	pDeath <- 1 - exp(-inds$Z*tincr)
	dead <- which(runif(nrow(inds)) < pDeath)
	# determine if natural or fished
	if(length(dead) > 0){
	  inds$alive[dead] <- 0
	  tmp <- cbind(inds$F[dead], inds$Z[dead])
	  # Fd=1 for fished individuals; Fd=0, for those that died naturally
	  Fd <- apply(tmp, 1, FUN=function(x){sample(c(0,1), 1, prob=c(M/x[2], x[1]/x[2]) )})
  	inds$Fd[dead] <- Fd
    rm(tmp)
	}
	return(inds)
}

remove.inds <- function(inds){
  dead <- which(inds$alive == 0)
  if(length(dead)>0) {inds <- inds[-dead,]}
  return(inds)
}

reproduce.inds <- function(inds){
	# reproduction can only occur of population contains >1 mature individual
	if(repro > 0 & sum(inds$mat) > 0){
		#calc. SSB
	  SSB <- sum(inds$W*inds$mat)
	  n.recruits <- ceiling(srrBH(rmax, beta, SSB) * repro)
		# make recruits 
		offspring <- make.inds(
			id = seq(lastID+1, length.out=n.recruits)
		)
		# express genes in recruits
		offspring <- express.inds(offspring)	
		#combine all individuals
		inds <- rbind(inds, offspring)
	}	
	return(inds)	
}

record.inds <- function(inds, ids=1:10, rec=NULL){
	if(is.null(rec)) {
		rec <- vector(mode="list", length(ids))
		names(rec) <- ids
		inds <- inds
	} else {
		ids <- as.numeric(names(rec))
	}
	if(length(rec) > 0) {
		inds.rows.rec <- which(!is.na(match(inds$id, ids)))
		if(length(inds.rows.rec) > 0){
			for(ii in inds.rows.rec){
				match.id <- match(inds$id[ii], ids)
				if(is.null(rec[[match.id]])) {
					rec[[match.id]] <- inds[ii,]
				} else {
					rec[[match.id]] <- rbind(rec[[match.id]], inds[ii,])
				}
			}
		}
	}
	rec
}




# run model ---------------------------------------------------------------

# Initial population
lastID <- 0
inds <- make.inds(
  id=seq(N0)
)
inds <- express.inds(inds)

# results object
res <- list()
res$pop <- list(
  dates = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) ),
	N = NaN*timeseq,
	B = NaN*timeseq,
	SSB = NaN*timeseq
)

# simulation
pb <- txtProgressBar(min=1, max=length(timeseq), style=3)
for(j in seq(timeseq)){
  t <- timeseq[j]
	
  # harvest rate applied? lfq sampled?
  if(length(fished_t) == 0){
    Fmax <- 0
    lfqSamp <- 0
  } else {
    if(min(sqrt((t-fished_t)^2)) < 1e-8){
      Fmax <- harvest_rate
      lfqSamp <- 1
    } else {
      Fmax <- 0
      lfqSamp <- 0
    }
  }
	
  repro <- repro_t[j]
	
	# population processes
	inds <- grow.inds(inds)
	inds <- mature.inds(inds)
	inds <- reproduce.inds(inds)
	inds <- death.inds(inds)
  if(lfqSamp){
    tmp <- try( sample(inds$L, ceiling(sum(inds$Fd)*lfqFrac), prob = inds$Fd), silent = TRUE)
    if(class(tmp) != "try-error"){lfq[[j]] <- tmp}
    rm(tmp) 
  }
	inds <- remove.inds(inds)
	
	# update results
	res$pop$N[j] <- nrow(inds)
	res$pop$B[j] <- sum(inds$W)
	res$pop$SSB[j] <- sum(inds$W*inds$mat)
	
	setTxtProgressBar(pb, j)

}
close(pb)



# Export data -------------------------------------------------------------

# Trim and Export 'lfq'
lfq2 <- lfq[which(sapply(lfq, length) > 0)]


# binned version of lfq
dates <- yeardec2date( date2yeardec(timemin.date) + (as.numeric(names(lfq2)) - timemin) )
Lran <- range(unlist(lfq2))
Lran[1] <- floor(Lran[1])
Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 1) * bin.size
bin.breaks <- seq(Lran[1], Lran[2], by=bin.size)
bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
res$lfqbin <- list(
  sample.no = seq(bin.mids),
  midLengths = bin.mids,
  dates = dates,
  catch = sapply(lfq2, FUN = function(x){
    hist(x, breaks=bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
  })
)

# record mean parameters
res$growthpars <- list(
  K = K.mu,
  Linf = Linf.mu,
  C = C,
  ts = ts,
  phiprime = phiprime,
  tmaxrecr = tmaxrecr
)

return(res)

} # end of function

