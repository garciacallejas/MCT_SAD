#' Calculate Hill's true diversity metric
#'
#' @param abundances vector of species abundances
#' @param q this parameter defines the kind of mean used (Tuomisto 2012)
#'
#' @return single diversity value
#' @export
#'
#' @examples hill.diversity(c(1,20,34,32,355,5))
hill.diversity <- function(abundances,q = 1){
  ab <- abundances[abundances != 0]
  if(length(ab)>0){
  ab <- ab/sum(ab)
  R <- length(ab)
  # hill diversity is not defined for q = 1,
  # but its limit exists and equals
  # the exponential of shannon entropy 
  # (Jost 2006,2007,Tuomisto 2012)
  if(q == 1){
    D <- exp(-sum(ab*log(ab)))
  }else{
    mean.p <- (sum(ab*ab^(q-1)))^(1/(q-1))
    D <- 1/mean.p
  }
  }else{
    D <- 0
  }
  return(D)
}