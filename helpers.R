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
  if(any(class(phenotypes) == 'list')) {
    all_names <- unique(unlist(purrr::map(phenotypes, ~colnames(.x))))
    index <- purrr::map(phenotypes, ~which(all_names %in% colnames(.x)))
    return(index)
  } else {
    index <- list(1:ncol(phenotypes))
    return(index)
  }
}
