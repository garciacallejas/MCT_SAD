
#' Ricker model with stronger density dependence: a_ii * Ni^2
#' with specific alpha values and no covariate effects
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra included for compatibility, not used in this model.
#' @param alpha_inter single numeric value.
#' @param lambda_cov included for compatibility, not used in this model.
#' @param alpha_cov included for compatibility, not used in this model.
#' @param abundance named numeric vector of abundances in the previous timestep.
#' @param covariates included for compatibility, not used in this model.
#'
#' @return numeric abundance projected one timestep
#' @export
HD_project_alpha_pairwise_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance,
                                                                   covariates){
  
  spnames <- names(abundance)
  
  # this should be ok
  # alpha <- c(alpha_intra,alpha_inter)
  # alpha <- alpha[spnames]
  
  numsp <- length(abundance)
  abund_intra <- abundance[which(names(abundance) == names(alpha_intra))]
  abund_inter <- abundance[which(names(abundance) %in% names(alpha_inter))]
  
  expected_abund <- NA_real_
  
  alpha_inter <- alpha_inter * -1
  alpha_intra <- alpha_intra * -1
  
  term <- 0
  for(z in 1:length(alpha_inter)){
    term <- term + alpha_inter[z]*abund_inter[z]
  }
  
  intra_term <- alpha_intra*abund_intra^2
  term <- term + intra_term
  
  # note the lambda inside the exponential
  # hopefully this prevents unbounded growth
  expected_abund <- (exp(lambda + term))*abund_intra
  
  expected_abund
}

