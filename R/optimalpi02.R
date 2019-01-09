#' my title
#'
#' @name optimalpi02
#'
#' @usage optimalpi02
#'
#' @keywords internal
#'
#' @return none

optimalpi02=function(n,ar,art,appv,aa,atau,lam1c=1.6,alam1e,lam0c=0.8,alam0e,precision=0.0001){
  value=10000
  pi0=0
  for (i in 0:(1/precision)){
    pi0tem=i*precision
    rstar=(art*appv+pi0tem*(ar-art*appv))/(art+(1-art)*pi0tem)
    ns=n/(art+(1-art)*pi0tem)
    p1c=1-aa/(ns*lam1c)*exp(-atau*lam1c)*(1-exp(-ns*lam1c/aa))
    p1e=1-aa/(ns*alam1e)*exp(-atau*alam1e)*(1-exp(-ns*alam1e/aa))
    p0c=1-aa/(ns*lam0c)*exp(-atau*lam0c)*(1-exp(-ns*lam0c/aa))
    p0e=1-aa/(ns*alam0e)*exp(-atau*alam0e)*(1-exp(-ns*alam0e/aa))
    v=2/(n*rstar)*(1/p1c+1/p1e)+2/(n*(1-rstar))*(1/p0c+1/p0e)
    if (v<value){
      value=v
      pi0=pi0tem
    }
  }
  return(pi0)
}
