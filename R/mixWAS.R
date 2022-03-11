#' Compute score vector and covariance matrix for a single site
#' @param snps: vector of SNPs in [0,1,2]
#' @param phenotypes: matrix of phenotypes (one per column), with names of phenotypes specified as column names
#' @param covariates: matrix or data frame of covariates
#' @param covariate_map: Default = NULL. If NULL, function assumes by all covariates are to be used for each phenotype.
#' If this is not desired behavior, user can supply a data frame with two columns
#' one called `variable` and a second called `phenotype`. In the variable column is the name of covariates, with phenotypes being specified as
#' 'all' (to use the variable for all phenotypes) or the name of a phenotype. Variables can be entered multiple times if they go to multiple phenotypes (but not all).
#' A phenotype specific data set of covariates with use 'all' covaraite + pheno-type sepcific covariates.
#' @param types: optional vector specifying data types ('continuous', 'binary', 'count'). Default = NULL (phenotype data types will be inferred)
#' Note that 'count' will never be inferred, only 'binary' or 'continuous'.
#' @return a list with components `score` and `V`, the score vector and covariance matrix for single site
#' @export
mixWAS_single_site <- function(snps, phenotypes, covariates, covariate_map = NULL, types = NULL) {
  ### Missing-ness (code as numeric to be 0 for each math but keep track of which are NA)
  Delta <- 1 - is.na(phenotypes)
  phenotypes[is.na(phenotypes)] <- 0

  ### Infer Phenotypes if types not supplied
  if(any(is.null(types))) {
    types <- unlist(infer_types(phenotypes))
  }

  ### number of phenotypes at site
  q <- ncol(phenotypes)

  ### Build list of all unique data sets needed
  ### Handle one hot encoding of character variables
  covariate_info <- compile_covariates(covariates, covariate_map, phenotypes, q)
  covariates <- covariate_info$covariates
  covariate_index <- covariate_info$index

  ### Run all regressions of reduced models
  gamma_list <- run_regressions(snps, phenotypes, covariates, covariate_index, Delta, types)

  ### Compute Score for all phenotypes at site
  score <- compute_site_score(snps, phenotypes, covariates, covariate_index, gamma_list, Delta, types)

  ### Compute computes necessary to build variance matrix
  ## Hessian
  hessian_list <- get_hessian_components(snps, covariates, covariate_index, gamma_list, types)

  ## Score Squared
  # First get each phenotype X/Z score components
  ss_components <- get_score_square_components(snps, phenotypes, covariates, covariate_index, gamma_list, Delta, types)

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


#' mixWAS Algorithm for All Sites (if all data can be supplied at once)
#' @param snps: list of snps (one for each site), each a vector of SNPs in [0,1,2]
#' @param phenotypes: list of phenotypes, each a matrix of phenotypes (one per column), with names of phenotypes specified as column names
#' @param covariates: list of covariates, each a matrix or data frame of covariates
#' @param covariate_map: default = NULL, list of covariate matrix build instructions for each site.
#' If specifying for all sites and just need some sites with complex build instruction (not all phenotypes use all covariates)
#' then supply a list of length = # of sites, entries = NULL for default behavior, and a data frame like object of instructions if needed.
#' See `mixWAS_single_site` help page for more information.
#' @param phenotype_index: list of vectors giving the index (numeric) of which phenotypes are in each site's matrix. If NULL (default), will be
#' infered from matrix colnames.
#' @param types: list of data types ('continuous', 'binary', 'count'). Default = NULL (infer types).
#' #' Note that 'count' will never be inferred, only 'binary' or 'continuous'.
#' @param parallel_sites: logical, if score/variance component computations should be parallelized over sites. Default = FALSE
#' @param return_p: logical, if TRUE return P-values, else return components like score/variance. Default = T
#' @return p-value of aggregate test
#' @export
mixWAS <- function(snps, phenotypes, covariates,
                   covariate_map = NULL, phenotype_index = NULL, types = NULL,
                   parallel_sites = F, return_p = T) {
  ### Handle the case of single site
  if(!is.list(snps)) {
    snps <- list(snps)
  }
  if(!is.list(phenotypes)) {
    phenotypes <- list(phenotypes)
  }
  if(!is.list(covariates)) {
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

  ### # of phenotypes
  q <- max(unlist(phenotype_index))

  ## Compute Score and Variance Matrix for Each Site
  if(parallel_sites) {
    mapper <- purrr::map
  } else {
    mapper <- furrr::future_map
  }

  testing_components <-
    mapper(1:num_sites, ~mixWAS_single_site(snps[[.x]],
                                            phenotypes[[.x]],
                                            covariates[[.x]],
                                            covariate_map[[.x]],
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
  z <- compute_z_scores(score, V_inv)

  ### Return components rather than p-values (helpful for parallelization)
  if(!return_p) {
    components <-
      list('score' = score,
           'V' = V,
           'V_inv' = V_inv,
           'z' = z,
           'q' = q)
    return(components)
  }

  ## Compute P-value(s)
  p_values <- run_hypothesis_tests(score, V_inv, z, q)

  return(p_values)
}

