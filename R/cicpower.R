#' my title
#'
#' @name cicpower
#'
#' @usage cicpower
#'
#' @keywords internal
#'
#' @return none

cicpower=function(rr,tref,n1,n2,m1c,m2c, r1,r2, accrual,tmin,alpha,pr){
  ## Find mortality in intervention groups
  if(m1c>1 | m2c>1)
    stop("m1c and m2c must be fractions")

  m1i <- (1-r1/100)*m1c
  m2i <- (1-r2/100)*m2c

  if(pr) {
    cat("\nAccrual duration:",accrual,"y  Minimum follow-up:",tmin,"y\n")
    cat("\nSample size Stratum 1:",n1,"  Stratum 2:",n2,"\n")
    cat("\nAlpha=",alpha,"\n")
    d <- list(c("Stratum 1","Stratum 2"), c("Control","Intervention"))
    m <- cbind(c(m1c,m2c),c(m1i,m2i))
    dimnames(m) <- d
    cat("\n",tref,"-year Mortalities\n",sep=""); print(m)
  }

  ## Find exponential hazards for all groups
  lam1c <- -logb(1-m1c)/tref
  lam2c <- -logb(1-m2c)/tref
  lam1i <- -logb(1-m1i)/tref
  lam2i <- -logb(1-m2i)/tref

  if(pr) {
    lam <- cbind(c(lam1c,lam2c),c(lam1i,lam2i))
    dimnames(lam) <- d
    cat("\nHazard Rates\n"); print(lam)
  }

  ## Find probability that a subject will have her event observed during
  ## the study, for all groups
  tmax <- tmin+accrual
  p1c <- 1-1/accrual/lam1c*(exp(-tmin*lam1c)-exp(-tmax*lam1c))
  p2c <- 1-1/accrual/lam2c*(exp(-tmin*lam2c)-exp(-tmax*lam2c))
  p1i <- 1-1/accrual/lam1i*(exp(-tmin*lam1i)-exp(-tmax*lam1i))
  p2i <- 1-1/accrual/lam2i*(exp(-tmin*lam2i)-exp(-tmax*lam2i))

  if(pr) {
    p <- cbind(c(p1c,p2c), c(p1i,p2i))
    dimnames(p) <- d
    cat("\nProbabilities of an Event During Study\n")
    print(p)
  }

  ##Find expected number of events, all groups
  m1c <- p1c*n1/2
  m2c <- p2c*n2/2
  m1i <- p1i*n1/2
  m2i <- p2i*n2/2

  if(pr) {
    m <- cbind(c(m1c,m2c), c(m1i,m2i))
    dimnames(m) <- d
    cat("\nExpected Number of Events\n")
    print(round(m,1))
  }

  ## Find expected value of observed log hazard ratio
  delta <- rr*logb(lam1i/lam1c)+(1-rr)*logb(lam2i/lam2c)
  if(pr)
    cat("\nRatio of hazard ratios:",format(exp(delta)),"\n")

  ## Find its variance
  v <- rr^2*(1/m1c + 1/m1i)+(1-rr)^2*(1/m2c + 1/m2i)
  sd <- sqrt(v)
  if(pr)
    cat("Standard deviation of log ratio of ratios:",format(sd),"\n")

  z <- -qnorm(alpha/2)
  ## if(pr) cat("\nCritical value:",format(z),"\n")

  c(Power = 1 - ( pnorm(z - abs(delta)/sd) - pnorm(-z - abs(delta)/sd) ) )
}
