
#' Ricker model for projecting abundances,
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
AP_project_alpha_pairwise_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance,
                                                                   covariates){
  
  
  # quickest way to include germination rates is to hard-code them here
  germ <- c(BEMA = 0.38, CETE = 0.53, CHFU = 0.8, CHMI = 0.76, HOMA = 0.94,
            LEMA = 0.89, LYTR = 0.62, MEEL = 0.59, MESU = 0.77, PAIN = 0.33,
            PLCO = 0.84, POMA = 0.85, POMO = 0.91, PUPA = 0.84, SASO = 0.52,
            SCLA = 0.69, SOAS = 0.6, SPRU = 0.44, SUSP = 0.61)
  
  surv <- c(BEMA = 0.43, CETE = 0.1, CHFU = 0.38, CHMI = 0.32, HOMA = 0.25,
            LEMA = 0.33, LYTR = 0.51, MEEL = 0.60, MESU = 0.63, PAIN = 0.60,
            PLCO = 0.57, POMA = 0.46, POMO = 0.49, PUPA = 0.55, SASO = 0.30,
            SCLA = 0.41, SOAS = 0.29, SPRU = 0.41, SUSP = 0.58)
  
  spnames <- names(abundance)
  
  focal.g <- germ[names(lambda)]
  focal.s <- surv[names(lambda)]
  
  # check germination rates against par names
  my.g <- germ[sort(spnames)]
  
  alpha <- c(alpha_intra,alpha_inter)
  alpha <- alpha[spnames]
  
  numsp <- length(abundance)
  expected_abund <- NA_real_
  
  num <- lambda
  term <- 0
  for(i.sp in 1:numsp){
    term <- term - alpha[i.sp]*log((my.g[i.sp]*abundance[i.sp]) + 1) 
  }# for each sp
  
  # term = 0 #create the denominator term for the model
  # for(z in 1:ncol(all_neigh_matrix)){
  #   term <- term - alpha[z]*log((my.g[z]*all_neigh_matrix[,z]) + 1) 
  # }
  # pred <- lambda * exp(term)
  
  expected_abund <- (lambda * exp(term)) * focal.g + (abundance[names(lambda)] * focal.s * (1-focal.g))
  expected_abund
}

