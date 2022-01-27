### Expit Function
expit <- function(x) {
  1/(1+exp(-x))
}

### Fast Diagonal  Matrix Multiplication
### When we have a diagonal matrix to multiply by another matrix we
### can do a lot faster than forming very large diagonal matrices
### d: vector of diagonal entries of matrix (A = diag(d)) but we never form this
### B: matrix to multiply by
### placement: 'left' or 'right', where the diagonal matrix A is placed relative to B
####           'left' = A %*% B
####           'right' = B %*% A
fast_diag_multiply <- function(d, B, placement = 'left') {
  if(placement == 'left') {
    M <- do.call(rbind, purrr::map(1:nrow(B), ~{d[.x] * B[.x,]}))
  } else {
    M <- do.call(cbind, purrr::map(1:ncol(B), ~{d[.x] * B[,.x]}))
  }
  return(M)
}

is_binary <- function(x) {
  x <- x[!is.na(x)]
  return(all(x == 0 | x == 1))
}


### Infer types for phenotype matrices for each site
###
### Will not ever infer count type. Will check for binary types otherwise
### will code phenotype as continuous
infer_types <- function(phenotypes) {
  if(any(class(phenotypes) == 'list')) {
    type_list <-
      purrr::map(phenotypes, ~ifelse(apply(.x, 2, is_binary), 'binary', 'continuous'))
    return(type_list)
  } else {
    types <- list(ifelse(apply(phenotypes, 2, is_binary), 'binary', 'continuous'))
    return(types)
  }
}

infer_phenotype_index <- function(phenotypes) {
  if(is.list(phenotypes)) {
    all_names <- unique(unlist(purrr::map(phenotypes, ~colnames(.x))))
    index <- purrr::map(phenotypes, ~which(all_names %in% colnames(.x)))
    return(index)
  } else {
    index <- list(1:ncol(phenotypes))
    return(index)
  }
}

### Functions to Build out list of covariate data sets for each phenotype
### as well as keeping track of which datasets correspond to which phenotype
compile_covariates <- function(covariates, covariate_map, phenotypes, q) {
  ### If there is no covariate map then the covariate list is just
  ### one data set for all
  if(is.null(covariate_map)) {
    index <- rep(1, q)
    covariates <- list(covariates)
  } else {
    unique_datasets <- unique(covariate_map$phenotype)

    ### List of varaibles in each data set
    variables <- purrr::map(unique_datasets, ~covariate_map$variable[covariate_map$phenotype == .x])
    if(length(variables) > 1) {
      variables <- c(variables[1], purrr::map(variables[-1], ~c(variables[[1]], .x)))
    }

    ### Build Data set
    covariates <- purrr::map(variables, ~covariates[,.x])
    index <- rep(1, q)
    ix <- colnames(phenotypes) %in% unique_datasets
    index[ix] <- purrr::map_dbl(colnames(phenotypes)[ix], ~which(unique_datasets == .x))
  }

  ### One Hot Encoding if Needed
  for(i in 1:length(covariates)) {
    data <- covariates[[i]]
    if(is.data.frame(data)) {
      covariates[[i]] <- as.matrix(fastDummies::dummy_cols(data, remove_first_dummy = T, remove_selected_columns = T))
    }
  }

  return(list('index' = index,
              'covariates' = covariates))
}
