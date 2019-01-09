#' my title
#'
#' @name SamSizes
#'
#' @usage SamSizes
#'
#' @keywords internal
#'
#' @return none

SamSizes=function(r,rt,ppv,a,total,tau,lam1c,lam1e,lam0c,lam0e,precision=0.0001,alpha,power0,N0=0){

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
    pow=Hmisc::ciapower(tref,N*r,N*(1-r),mc[1],mc[2],red[1],red[2],accrual,total-accrual,
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
    N=N+10
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
  follow1=r/2*min(max(0,total-1/lam1c),N1/a)*a*1/lam1c+max(0,N1/a-max(0,total-1/lam1c))*a*r/2*(1/2*max(0,N1/a-max(0,total-1/lam1c))+total-N1/a)
  follow2=r/2*min(max(0,total-1/lam1e),N1/a)*a*1/lam1e+max(0,N1/a-max(0,total-1/lam1e))*a*r/2*(1/2*max(0,N1/a-max(0,total-1/lam1e))+total-N1/a)
  follow3=(1-r)/2*min(max(0,total-1/lam0c),N1/a)*a*1/lam0c+max(0,N1/a-max(0,total-1/lam0c))*a*(1-r)/2*(1/2*max(0,N1/a-max(0,total-1/lam0c))+total-N1/a)
  follow4=(1-r)/2*min(max(0,total-1/lam0e),N1/a)*a*1/lam0e+max(0,N1/a-max(0,total-1/lam0e))*a*(1-r)/2*(1/2*max(0,N1/a-max(0,total-1/lam0e))+total-N1/a)
  follow=follow1+follow2+follow3+follow4
  cost1=5000*N1+200*follow+1000*N1;
  rs=(rt*ppv+pi0*(r-rt*ppv))/(rt+(1-rt)*pi0)
  ns=ceiling(ns)
  follow1=rs/2*min(max(0,total-1/lam1c),ns/a)*a*1/lam1c+max(0,ns/a-max(0,total-1/lam1c))*a*rs/2*(1/2*max(0,ns/a-max(0,total-1/lam1c))+total-ns/a)
  follow2=rs/2*min(max(0,total-1/lam1e),ns/a)*a*1/lam1e+max(0,ns/a-max(0,total-1/lam1e))*a*rs/2*(1/2*max(0,ns/a-max(0,total-1/lam1e))+total-ns/a)
  follow3=(1-rs)/2*min(max(0,total-1/lam0c),ns/a)*a*1/lam0c+max(0,ns/a-max(0,total-1/lam0c))*a*(1-rs)/2*(1/2*max(0,ns/a-max(0,total-1/lam0c))+total-ns/a)
  follow4=(1-rs)/2*min(max(0,total-1/lam0e),ns/a)*a*1/lam0e+max(0,ns/a-max(0,total-1/lam0e))*a*(1-rs)/2*(1/2*max(0,ns/a-max(0,total-1/lam0e))+total-ns/a)
  follow=follow1+follow2+follow3+follow4
  cost2=100*ns+5000*N2+200*follow+1000*N2;

  return(c(N1,N2,ns,pi0,cost1,cost2))
}
