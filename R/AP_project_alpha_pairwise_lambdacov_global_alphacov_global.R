
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
AP_project_alpha_pairwise_lambdacov_global_alphacov_global <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance,
                                                                   covariates){
  
  
  ##################
  # lambda <- c(BEMA = 10.00102)
  # alpha_intra <- c(LEMA = -0.0009885923)
  # alpha_inter <- c(BEMA = -0.0009948755,CHFU = -0.0009844057,HOMA = 0.0001998094)
  # lambda_cov <- 0.0195262
  # alpha_cov <- list(precipitation = 0)
  # abundance <- c(BEMA = 2,CHFU = 14,HOMA = 44, LEMA = 10)
  # covariates <- 625.5
  
  
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
  # alpha_cov_order <- substr(names(alpha_cov),nchar(names(alpha_cov))-3,nchar(names(alpha_cov)))
  # alpha_cov <- alpha_cov[match(names(my.g),alpha_cov_order)]
  #
  
  # model
  num = 1
  focal.cov.matrix <- as.matrix(covariates)
  for(v in 1:ncol(focal.cov.matrix)){
    num <- num + lambda_cov[v]*focal.cov.matrix[,v] 
  }
  
  # denominator
  cov_term <- 0 
  for(v in 1:ncol(focal.cov.matrix)){
    cov_term <- cov_term + alpha_cov[[v]] * focal.cov.matrix[,v]
  }
  
  numsp <- length(abundance)
  expected_abund <- NA_real_
  
  term <- 0 #create the denominator term for the model
  for(z in 1:length(abundance)){
    term <- term - (alpha[z] + cov_term)*log((my.g[z]*abundance[z]) + 1) 
  }
  expected_abund <- (lambda * (num) * exp(term)) * focal.g + (abundance[names(lambda)] * focal.s * (1-focal.g))
  expected_abund
  
  # num <- lambda
  # term <- 0
  # for(i.sp in 1:numsp){
  #   term <- term - alpha[i.sp]*log((my.g[i.sp]*abundance[i.sp]) + 1) 
  # }# for each sp
  # 
  # expected_abund <- (lambda * exp(term)) * focal.g + (abundance[names(lambda)] * focal.s * (1-focal.g))
  # expected_abund
  
}

