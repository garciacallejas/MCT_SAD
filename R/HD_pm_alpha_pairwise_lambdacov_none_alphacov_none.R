#' Ricker model with stronger density dependence: a_ii * Ni^2
#'
#' @param par 1d vector of initial parameters: 'lambda', 'alpha_intra' (optional), 'alpha_inter', and 'sigma'
#' @param fitness 1d vector of fitness observations, in log scale
#' @param neigh_intra_matrix optional matrix of one column, number of intraspecific neighbours for each observation
#' @param neigh_inter_matrix matrix of arbitrary columns, number of interspecific neighbours for each observation
#' @param covariates included for compatibility, not used in this model
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","alpha_intra","alpha_inter".
#'
#' @return log-likelihood value
#' @export
HD_pm_alpha_pairwise_lambdacov_none_alphacov_none <- function(par,
                                                              fitness,
                                                              neigh_intra_matrix = NULL,
                                                              neigh_inter_matrix,
                                                              covariates,
                                                              fixed_parameters){

  
  # retrieve parameters -----------------------------------------------------
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha_intra,alpha_inter,alpha_cov,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  if(is.null(fixed_parameters[["lambda"]])){
    lambda <- par[pos]
    pos <- pos + 1
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # if(is.null(fixed_parameters$lambda_cov)){
  #   lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
  #   pos <- pos + ncol(covariates)
  # }else{
  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
  # }
  
  if(!is.null(neigh_intra_matrix)){
    # intra
    if(is.null(fixed_parameters[["alpha_intra"]])){
      alpha_intra <- par[pos]
      pos <- pos + 1
    }else{
      alpha_intra <- fixed_parameters[["alpha_intra"]]
    }
  }else{
    alpha_intra <- NULL
  }
  
  # inter
  if(is.null(fixed_parameters[["alpha_inter"]])){
    alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
    pos <- pos + ncol(neigh_inter_matrix) -1
  }else{
    alpha_inter <- fixed_parameters[["alpha_inter"]]
  }
  
  # if(is.null(fixed_parameters$alpha_cov)){
  #   alpha.cov <- par[pos:(pos+(ncol(covariates)*ncol(neigh_matrix))-1)]
  #   pos <- pos + (ncol(covariates)*ncol(neigh_matrix))
  # }else{
  #   alpha.cov <- fixed_parameters[["alpha.cov"]]
  # }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  # if(!is.null(alpha_intra)){
  #   alpha <- c(alpha_intra,alpha_inter)
  #   all_neigh_matrix <- cbind(neigh_intra_matrix,neigh_inter_matrix)
  # }else{
  #   alpha <- alpha_inter
  #   all_neigh_matrix <- neigh_inter_matrix
  # }
  
  alpha_inter <- alpha_inter * -1
  alpha_intra <- alpha_intra * -1
  
  term = 0 #create the denominator term for the model
  for(z in 1:ncol(neigh_inter_matrix)){
    term <- term + alpha_inter[z]*neigh_inter_matrix[,z] 
  }
  
  # neigh_intra_matrix is always of dimensions [n,1]
  # so this keeps it as a vector
  intra_term <- alpha_intra*neigh_intra_matrix[,1]^2
  term <- term + intra_term
  
  # note the lambda inside the exponential
  # hopefully this prevents unbounded growth
  pred <- exp(lambda + term)

  if(any(pred == 0)){
    pred[which(pred == 0)] <- 1e-10
  }
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}