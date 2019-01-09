#' my title
#'
#' @name recalcuint
#'
#' @usage recalcuint
#'
#' @keywords internal
#'
#' @return none

recalcuint=function(N2A,nsA,N2B,nsB,pi0B,ppv,r,rt){
  rsB=(rt*ppv+pi0B*(r-rt*ppv))/(rt+(1-rt)*pi0B);
  if (N2A*ppv<=N2B*rsB){
    N=N2B;rs=rsB;ns=nsB;pi0=pi0B;
  } else {
    r2=N2A*ppv/(N2B*(1-rsB));
    if (r2>=ppv/(1-ppv)) {
      N=N2A; rs=ppv;
      ns=nsA; pi0=0;
    } else {
      pi0=(rt*ppv-r2*rt*(1-ppv))/(r2*(1-rt-r+rt*ppv)-r+rt*ppv);
      rs=(rt*ppv+pi0*(r-rt*ppv))/(rt+(1-rt)*pi0);
      N=N2A*ppv/rs;
      ns=N/(rt+(1-rt)*pi0);
    }
  }
  return(c(N,rs,ns,pi0))
}
