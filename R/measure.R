#' my title
#'
#' @name measure
#'
#' @usage measure
#'
#' @keywords internal
#'
#' @return none

measure=function(lam1e,lam0e,lam1c,lam0c,total,r,N1,a,N2,rs,ns,rt,ppv,costscr,costtrt,costfol,costtest){
  follow1=r/2*min(max(0,total-1/lam1c),N1/a)*a*1/lam1c+max(0,N1/a-max(0,total-1/lam1c))*a*r/2*(1/2*max(0,N1/a-max(0,total-1/lam1c))+total-N1/a)
  follow2=r/2*min(max(0,total-1/lam1e),N1/a)*a*1/lam1e+max(0,N1/a-max(0,total-1/lam1e))*a*r/2*(1/2*max(0,N1/a-max(0,total-1/lam1e))+total-N1/a)
  follow3=(1-r)/2*min(max(0,total-1/lam0c),N1/a)*a*1/lam0c+max(0,N1/a-max(0,total-1/lam0c))*a*(1-r)/2*(1/2*max(0,N1/a-max(0,total-1/lam0c))+total-N1/a)
  follow4=(1-r)/2*min(max(0,total-1/lam0e),N1/a)*a*1/lam0e+max(0,N1/a-max(0,total-1/lam0e))*a*(1-r)/2*(1/2*max(0,N1/a-max(0,total-1/lam0e))+total-N1/a)
  follow=follow1+follow2+follow3+follow4
  cost1=costtrt*N1+costfol*follow+costtest*N1;
  follow1=rs/2*min(max(0,total-1/lam1c),ns/a)*a*1/lam1c+max(0,ns/a-max(0,total-1/lam1c))*a*rs/2*(1/2*max(0,ns/a-max(0,total-1/lam1c))+total-ns/a)
  follow2=rs/2*min(max(0,total-1/lam1e),ns/a)*a*1/lam1e+max(0,ns/a-max(0,total-1/lam1e))*a*rs/2*(1/2*max(0,ns/a-max(0,total-1/lam1e))+total-ns/a)
  follow3=(1-rs)/2*min(max(0,total-1/lam0c),ns/a)*a*1/lam0c+max(0,ns/a-max(0,total-1/lam0c))*a*(1-rs)/2*(1/2*max(0,ns/a-max(0,total-1/lam0c))+total-ns/a)
  follow4=(1-rs)/2*min(max(0,total-1/lam0e),ns/a)*a*1/lam0e+max(0,ns/a-max(0,total-1/lam0e))*a*(1-rs)/2*(1/2*max(0,ns/a-max(0,total-1/lam0e))+total-ns/a)
  follow=follow1+follow2+follow3+follow4
  cost2=costscr*ns+costtrt*N2+costfol*follow+costtest*N2;
  pi00=(rt*ppv-rs*rt)/(rs-rs*rt-r+rt*ppv)
  return(c(pi00,N2,ns,N1,N2/N1,cost2/cost1,ns/N1))
}
