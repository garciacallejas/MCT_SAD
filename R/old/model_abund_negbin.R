#' Project abundance of individuals according to the Beverton-Holt third model
#'
#' @param sp.par dataframe with species in rows, and the following columns:
#' lambda: fecundity term
#' germ.rate: seed germination rate
#' survival.rate: annual survival of ungerminated seed
#' @param init.abund number of individuals at time t
#' @param cov.values Not used in model_abundBH3
#' @param alpha.matrix competition matrix
#' @param lambda.cov.matrix Not used in model_abundBH3
#' @param alpha.cov.matrix Not used in model_abundBH3
#' @param return.seeds boolean flag, whether the prediction should return 
#' number of seeds (i.e. \eqn{N_{i,t+1}}, eq. 1 of Lanuza et al. 2018), or number of
#' adult individuals, (i.e. \eqn{N_{i,t+1} * g} )
#'
#' @return 1d vector with number of individuals of each species at time t+1
#' @export
model_abund_negbin <- function(sp.par,init.abund,cov.values,alpha.matrix,lambda.cov.matrix,alpha.cov.matrix,return.seeds = FALSE){
  expected.abund <- rep(0,nrow(sp.par))
  for(i.sp in 1:nrow(sp.par)){
    
    interactions <- 0
    
    # TODO
    # for now, only adults in init.abund are considered
    for(j.sp in 1:nrow(sp.par)){
      interactions <- interactions + alpha.matrix[i.sp,j.sp]*init.abund[j.sp]
    }# for j.sp
    
    exponent <- sp.par$lambda[i.sp] - interactions
    fitness <- exp(exponent)
    
    expected.abund[i.sp] <- ((1-sp.par$germ.rate[i.sp])*sp.par$survival.rate[i.sp]) + sp.par$germ.rate[i.sp]*fitness
  }
  # just in case
  expected.abund[expected.abund < 0] <- 0
  # the above formulation gives number of seeds. If necessary, multiply by the germination rate to get number of individuals
  if(return.seeds){
    return(expected.abund)
  }else{
    return(expected.abund*sp.par$germ.rate)
  }
}

################

# similar function but projecting abundances, instead of returning -loglik

