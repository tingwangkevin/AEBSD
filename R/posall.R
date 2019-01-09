#' Sample Size of AEBSD for Postive and Overall
#'
#' Calculate minimum sample size needed of AEBSD for testing both treatment effect in biomarker positve and
#' overall treatment effect, and its comparison with BSD
#'
#' @name posall
#'
#' @usage posall(rtrue,raux,ppv,lam1e,lam0e,lam1c,lam0c,tau,raccrual,powerdesired,alpha,costs)
#'
#' @param rtrue the prevalence rates of patients with positive true biomarker
#' @param raux  the prevalence rates of patients with positive auxiliary biomarker
#' @param lam1e hazrd rate for biomarker-positive subpopulation in the experimental treatment group
#' @param lam0e hazrd rate for biomarker-negative subpopulation in the experimental treatment group
#' @param lam1c hazrd rate for biomarker-positive subpopulation in the control treatment group
#' @param lam0c hazrd rate for biomarker-negative subpopulation in the control treatment group
#' @param tau   minimum follow-up time
#' @param raccrual accrual rate
#' @param powerdesired desired design power
#' @param alpha    type-I error (for each test)
#' @param costs cost vector including costscr (scrrening cost), costtrt (average treatment cost), costfol
#' (follow-up cost) and costtest (biomarker test cost) for each randomized patient
#'
#' @return (1) The probability of randomizing patients with positive and negative auxiliary variable;
#' (2) Number of randomized patients for AEBSD (3) Screening number for AEBSD
#' (4) Number of randomized patients for BSD (5) Patient ratio (6) Cost ratio (7) Screening ratio
#'
#' @export

options(scipen=200)
posall=function(rtrue,raux,ppv,lam1e,lam0e,lam1c,lam0c,tau,raccrual,powerdesired,alpha,costs){
  costscr=costs[1]; costtrt=costs[2]; costfol=costs[3]; costtest=costs[4]
  r=rtrue; rt=raux; a=raccrual; power0=powerdesired
  lamc=c(lam1c,lam0c)
  lame=c(lam1e,lam0e)
  N2=10^100; pi0o=-0.001
  for (j in 1:1001){
    print(j-1)
    pi0=(j-1)/1000
    #cib
    pow=0; N=1
    red=c();mc=c();me=c()
    while(pow<power0){
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
      pow=cibpower(tref,n1,n2,mc[1],mc[2],red[1],red[2],accrual,tau,alpha,pr=FALSE)
      N=N+1
    }
    Nb=N
    #Nb=0
    #cic
    pow=0;N=1
    red=c();mc=c();me=c()
    while(pow<power0){
      #pi0=optimalpi02(n=N,ar=r,art=rt,appv=ppv,alam1e=lam1e,alam0e=lam0e,atau=tau,aa=a)
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
      pow=cicpower(r,tref,n1,n2,mc[1],mc[2],red[1],red[2],accrual,tau,alpha,pr=FALSE)
      N=N+1
    }
    Na=N
    N2new=max(Na,Nb)
    if(N2new<N2){
      pi0o=pi0
      N2=N2new
      Nbmin=Nb
      Namin=Na
    }
  }

  ns=N2/(rt+(1-rt)*pi0o)
  rs=(rt*ppv+pi0o*(r-rt*ppv))/(rt+(1-rt)*pi0o);
  total=ns/a+tau;##By the Ns we get, determine the total time in which follow up time of BED is 1-2 years
  C=positive1(ppv=ppv,lam1e=lam1e,lam0e=lam0e,lam1c=lam1c,lam0c=lam0c,r=r,rt=rt,total=total,tau=tau,a=a,alpha = 0.025,power0=0.9);##Get the new sample sizes for two methods
  pow=0;N=1
  red=c();mc=c();me=c()
  while(pow<power0){
    accrual=N/a
    tref=1
    for(i in 1:2){
      mc[i]=1-exp(-tref*lamc[i])
      me[i]=1-exp(-tref*lame[i])
      red[i]=(mc[i]-me[i])/mc[i]*100
    }
    n1=N*r
    n2=N*(1-r)
    pow=cicpower(r,tref,n1,n2,mc[1],mc[2],red[1],red[2],accrual,total-accrual,alpha,pr=FALSE)
    N=N+1
  }
  N1c=N
  N1=max(C[1],N1c)
  if (Namin>C[1]){
    D=positive0(ppv=ppv,lam1e=lam1e,lam0e=lam0e,lam1c=lam1c,lam0c=lam0c,r=r,rt=rt,tau=tau,a=a,alpha = 0.025,power0=0.9);##Get the new sample sizes for two methods
    pow=0;N=D[1]
    red=c();mc=c();me=c()
    while(pow<power0){
      accrual=N/a
      tref=1
      for(i in 1:2){
        mc[i]=1-exp(-tref*lamc[i])
        me[i]=1-exp(-tref*lame[i])
        red[i]=(mc[i]-me[i])/mc[i]*100
      }
      n1=N*r
      n2=N*(1-r)
      pow=cicpower(r,tref,n1,n2,mc[1],mc[2],red[1],red[2],accrual,tau,alpha,pr=FALSE)
      N=N+1
    }
    N1=N
    N2=N
    rs=r
    ns=N2
    total=N1/a+tau
  }

  E=measure(lam1e=lam1e,lam0e=lam0e,lam1c=lam1c,lam0c=lam0c,total=total,r=r,N1=N1,a=a,N2=N2,rs=rs,ns=ns,rt=rt,ppv=ppv);
  #E=c(Nbmin,Namin,C[1],N1c,N2,ns,N1)
  c(E) ##Each row contains (PPV,r,rt,N_BSD,N_BED,Ns,pi0,Cost_BSD,Cost_BED)
}
