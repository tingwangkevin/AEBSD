#' my title
#'
#' @name SamSizeso
#'
#' @usage SamSizeso
#'
#' @keywords internal
#'
#' @return none

SamSizeso=function(r,rt,ppv,a,tau,lam1c,lam1e,lam0c,lam0e,precision=0.0001,alpha,power0,N0=0){

  ################################################
  #numerically find the sample size satisfies######
  #centain power######

  #r:percentage of positive marker
  #rt:percentage of positive surrogate marker
  #ppv:percentage of positive marker when idt is positive surrogate marker
  #a:accrual rate (numbers of patients recruited) in yr
  #tau::min time to follow up in yr
  #precision:precision of pi0
  #power0:target power
  #pow:current power
  #N:total number of subjects for testing marker(1:1)
  #ns: screen number for testing surrogate marker
  #red:%reduction
  #tref:time point for comparison

  lamc=c(lam1c,lam0c)
  lame=c(lam1e,lam0e)

  #BSD
  pow=0
  N=N0

  red=c()
  mc=c()
  me=c()

  while(pow<power0){
    N=N+10
    accrual=N/a
    tref=1
    for(i in 1:2){
      mc[i]=1-exp(-tref*lamc[i])
      me[i]=1-exp(-tref*lame[i])
      red[i]=(mc[i]-me[i])/mc[i]*100
    }
    pow=Hmisc::ciapower(tref,N*r,N*(1-r),mc[1],mc[2],red[1],red[2],accrual,tau,
                 alpha,pr=FALSE)
  }
  N1=N

  #BED
  pow=0
  N=N0

  red=c()
  mc=c()
  me=c()

  while(pow<power0){
    N=N+2
    pi0=optimalpi02(n=N,ar=r,art=rt,appv=ppv,alam1e=lam1e,alam0e=lam0e,atau=tau,aa=a)
    ns=N/(rt+(1-rt)*pi0)
    accrual=ns/a
    tref=1
    for(i in 1:2){
      mc[i]=1-exp(-tref*lamc[i])
      me[i]=1-exp(-tref*lame[i])
      red[i]=(mc[i]-me[i])/mc[i]*100
    }
    n1=ns*(rt*ppv+pi0*(r-rt*ppv))
    n2=N-n1
    pow=Hmisc::ciapower(tref,n1,n2,mc[1],mc[2],red[1],red[2],accrual,tau,
                 alpha,pr=FALSE)
  }
  N2=N

  return(c(N1,N2,ns,pi0))
}
