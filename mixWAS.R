source('helpers.R')
source('run_regressions.R')
source('variance_components.R')
source('hypothesis_tests.R')

### Compute Variance and Score Component for single site
###
### snps: vector of SNPs in [0,1,2]
### phenotypes: Matrix of phenotypes (one per column)
### covariates: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates
### types: vector specifies the type of phenotype (e.g. continuous, binary, count), will be infered if NULL
mixWAS_single_site <- function(snps, phenotypes, covariates, types = NULL) {
  ### Missing-ness (code as numeric to be 0 for each math but keep track of which are NA)
  Delta <- 1 - is.na(phenotypes) 
  phenotypes[is.na(phenotypes)] <- 0
  
  ### Infer Phenotypes if types not supplied
  if(any(is.null(types))) {
    types <- infer_types(phenotypes)
  }
  
  ### number of phenotypes at site
  q <- ncol(phenotypes) 
  
  ### Run all regressions of reduced models
  gamma_list <- run_regressions(snps, phenotypes, covariates, Delta, types)
  
  ### Compute Score for all phenotypes at site
  score <- compute_site_score(snps, phenotypes, covariates, gamma_list, Delta, types)
  
  ### Compute computes necessary to build variance matrix
  ## Hessian
  hessian_list <- get_hessian_components(snps, covariates, gamma_list, types)
  
  ## Score Squared
  # First get each phenotype X/Z score components
  ss_components <- get_score_square_components(snps, phenotypes, covariates, gamma_list, Delta, types)
  
  # Next actually compute the score squared for each unique pair
  indices <- 
    tidyr::crossing('i' = 1:q, 'j' = 1:q) %>% 
    dplyr::filter(i <= j) %>% 
    dplyr::mutate('row_id' = 1:nrow(.))
  
  ss_list <- build_score_square_list(indices, ss_components)
  
  
  ### Build Variance Matrix
  V <- build_variance_matrix(ss_list, hessian_list, indices, q)
  
  return(list('score' = score, 'V' = V))
}


### mixWAS algoirthm
###
### snps: vector of SNPs in [0,1,2]
### phenotypes: list of matrices (one per site) of phenotypes (one phenotype per column)
### covariates: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates
### types: vector specifies the type of phenotype (e.g. continuous, binary, count), will be inferred if NULL
mixWAS <- function(snps, phenotypes, covariates, phenotype_index = NULL, types = NULL) {
  ### Handle the case of single site 
  if(all(class(snps) != 'list')) {
    snps <- list(snps)
  } 
  if(all(class(phenotypes) != 'list')) {
    phenotypes <- list(phenotypes)
  }
  if(all(class(covariates) != 'list')) {
    covariates <- list(covariates)
  }
  
  ### # of sites
  num_sites <- length(snps)
  
  ### Infer Phenotypes if types not supplied
  if(any(is.null(types))) {
    types <- infer_types(phenotypes)
  }
  
  ### Infer Phenotype index from names if not provided
  if(any(is.null(phenotype_index))) {
    phenotype_index <- infer_phenotype_index(phenotypes)
  }
  
  
  ### Compute Score and Variance Matrix for Each Site
  testing_components <- 
    purrr::map(1:num_sites, ~mixWAS_single_site(snps[[.x]],
                                                phenotypes[[.x]],
                                                covariates[[.x]],
                                                types[[.x]]))
  
  ### Build combined score/variance matrix
  score <- rep(0, q)
  V <- matrix(0, nrow = q, ncol = q)
  
  for(i in 1:num_sites) {
    index <- phenotype_index[[i]] 
    score[index] <- score[index] + testing_components[[i]]$score
    V[index, index] <- V[index, index] + testing_components[[i]]$V
  }
  
  V_inv <- solve(V)
  
  ### Compute P-value for SNP
  p1 <- score_test(score, V_inv, q)
  p2 <- max_text(score, V_inv, q)
  p_snp <- acat(c(p1, p2))
  
  return(p_snp)
}