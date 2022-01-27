### Functions to compute the various components 
### necessary to build variance matrix

### Compute Score for a Single Phenotype
###
### x: vector of SNPs in [0,1,2]
### y: vector of phenotypes 
### Z: matrix of site-specific, pheno-type specific covariates
### gamma: vector of coefficients for the Z covariates
### delta: vector indicating if phenotype is present (1) or missing (0)
### type: specifies the type of phenotype (e.g. continuous, binary, count)

compute_phenotype_score <- function(x, y, Z, gamma, delta, type) {
  ### Augment Covariates w/ intercept
  design <- cbind(1, Z)
  
  if(type == 'continuous') {
    score <- sum(x * (y - (design %*% gamma)) * delta)
  } else if(type == 'binary') {
    score <- sum(x * (y - expit(design %*% gamma)) * delta)
  } else if(type == 'count') {
    score <- sum(x * (y - exp(design %*% gamma)) * delta)
  }
  
  return(score)
}

### Compute Score for all Phenotypes with a single site
###
### x: vector of SNPs in [0,1,2]
### Y: Matrix of phenotypes (one per column)
### Z_list: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates
### gamma_list: list of coefficient vectors, each corresponding to a set of site-specific, pheno-type specific covariates
### Delta: vector indicating if phenotype is present (1) or missing (0)
### types: vector specifies the type of phenotype (e.g. continuous, binary, count)

compute_site_score <- function(x, Y, Z_list, Z_index, gamma_list, Delta, types) {
  score_vec <- 
    purrr::map_dbl(1:ncol(Y), ~compute_phenotype_score(x = x, 
                                                       y = Y[,.x],
                                                       Z = Z_list[[Z_index[.x]]],
                                                       gamma = gamma_list[[.x]],
                                                       delta = Delta[,.x],
                                                       type = types[.x])) 
  
  return(score_vec)
}

### Compute X/Z Score Vectors for a Single Phenotype 
##  (this will be used to compute score squared)
###
### x: vector of SNPs in [0,1,2]
### y: vector of phenotypes 
### Z: matrix of site-specific, pheno-type specific covariates
### gamma: vector of coefficients for the Z covariates
### delta: vector indicating if phenotype is present (1) or missing (0)
### type: specifies the type of phenotype (e.g. continuous, binary, count)

score_squared_components <- function(x, y, Z, gamma, delta, type) {
  design <- cbind(1, Z)
  
  ### Score Components for Phenotype 
  if(type == 'continuous') {
    score_x <- x * (y - (design %*% gamma)) * delta
    score_z <- fast_diag_multiply(d = as.vector((y-(design %*% gamma)) * delta),
                                  B = design,
                                  placement = 'left')
  } else if(type == 'binary') {
    score_x <- x * (y - expit(design %*% gamma)) * delta
    score_z <- fast_diag_multiply(d = as.vector((y-expit(design %*% gamma)) * delta),
                                  B = design,
                                  placement = 'left')
  } else if(type == 'count') {
    score_x <- x * (y - exp(design %*% gamma)) * delta
    score_z <- fast_diag_multiply(d = as.vector((y-exp(design %*% gamma)) * delta),
                                  B = design,
                                  placement = 'left')
  }
  
  score <- list('x' = score_x,
                'z' = score_z)
  
  return(score)
}

### Compute Score Squared Building Blocks 
### 
### x: vector of SNPs in [0,1,2]
### Y: Matrix of phenotypes (one per column)
### Z_list: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates
### gamma_list: list of coefficient vectors, each corresponding to a set of site-specific, pheno-type specific covariates
### Delta: Matrix vector indicating if phenotype is present (1) or missing (0)
### types: vector specifies the type of phenotype (e.g. continuous, binary, count)

get_score_square_components <- function(x, Y, Z_list, Z_index, gamma_list, Delta, types) {
  ss_components <- 
    purrr::map(1:ncol(Y), ~score_squared_components(x = x,
                                                    y = Y[,.x],
                                                    Z = Z_list[[Z_index[.x]]],
                                                    gamma = gamma_list[[.x]],
                                                    delta = Delta[,.x],
                                                    type = types[.x]))
  
  return(ss_components)
  
}

### Carry out the multiplication for the squared score component
square_score <- function(score_list1, score_list2) {
  ### Unpack score lists
  score_x1 <- score_list1$x
  score_z1 <- score_list1$z
  score_x2 <- score_list2$x
  score_z2 <- score_list2$z
  
  ### Score Square terms
  score_square_xx <- sum(score_x1 * score_x2)
  score_square_xz1 <- t(score_x1) %*% score_z2
  score_square_xz2 <- t(score_x2) %*% score_z1
  score_square_zz <- t(score_z1) %*% score_z2
  
  score_squared <- list('xx' = score_square_xx,
                        'xz1' = score_square_xz1,
                        'xz2' = score_square_xz2,
                        'zz' = score_square_zz)
  
  return(score_squared)
}

### Function to carry out all of the multiplication for the squared score component
build_score_square_list <- function(indices, ss_components) {
  ss_list <- 
    purrr::map2(indices$i, 
                indices$j,
                ~square_score(ss_components[[.x]], ss_components[[.y]])) 
  return(ss_list)
}

### x: vector of SNPs in [0,1,2]
### Z: matrix of site-specific, pheno-type specific covariates
### gamma: vector of coefficients for the Z covariates
### type: specifies the type of phenotype (e.g. continuous, binary, count)
compute_hessian <- function(x, Z, gamma, type) {
  design <- cbind(1, Z)
  
  if(type == 'continuous') {
    hessian_xz <-  t(x) %*% design
    hessian_zz <- t(design) %*% design
  } else if(type == 'binary') {
    ### pre-compute the left side of equation using our fast diagonal 
    ### matrix multiplication
    d <- as.vector((1-expit(design %*% gamma)) * expit(design %*% gamma))
    left_side <- fast_diag_multiply(d, design, 'left') 
    
    hessian_zz <- t(design) %*% left_side
    hessian_xz <- t(x) %*% left_side
  } else if(type == 'count') {
    ### pre-compute the left side of equation using our fast diagonal 
    ### matrix multiplication
    d <- as.vector(exp(design %*% gamma))
    left_side <- fast_diag_multiply(d, design, 'left') 
    
    hessian_zz <- t(design) %*% left_side
    hessian_xz <- t(x) %*% left_side
  }
  
  hessian <- list('xz' = hessian_xz,
                  'zz' = hessian_zz,
                  'zz_inv' = solve(hessian_zz))
  
  return(hessian)
}

### Compute Score for all Phenotypes with a single site
###
### x: vector of SNPs in [0,1,2]
### Z_list: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates
### gamma_list: list of coefficient vectors, each corresponding to a set of site-specific, pheno-type specific covariates
### types: vector specifies the type of phenotype (e.g. continuous, binary, count)

get_hessian_components <- function(x, Z_list, Z_index, gamma_list, types) {
  hessian_list <- 
    purrr::map(1:length(gamma_list), ~compute_hessian(x = x, 
                                                  Z = Z_list[[Z_index[.x]]],
                                                  gamma = gamma_list[[.x]],
                                                  type = types[.x])) 
  return(hessian_list)
  
}

### Compute the variance matrix entry v_ij given the pre-computed components
### for the score_squared for phenotypes ixj, hessian for phenotype i and phenotype j
compute_variance_entry <- function(score_squared, hessian_1, hessian_2) {
  ### Unpack pre-computed components
  h_xz_1 <- hessian_1$xz
  h_zz_inv_1 <- hessian_1$zz_inv
  h_xz_2 <- hessian_2$xz
  h_zz_inv_2 <- hessian_2$zz_inv
  
  ss_xx <- score_squared$xx
  ss_xz_1 <- score_squared$xz1
  ss_xz_2 <- score_squared$xz2
  ss_zz <- score_squared$zz
  
  ### Compute component
  v_ij <- 
    ss_xx - 
    h_xz_1 %*% h_zz_inv_1 %*% t(ss_xz_2) -  
    h_xz_2 %*% h_zz_inv_2 %*% t(ss_xz_1) + 
    h_xz_1 %*%  h_zz_inv_1 %*% ss_zz %*% h_zz_inv_2 %*% t(h_xz_2)
  
  return(v_ij)
}

### Build Variance Matrix
###
### ss_list: list of all pre-computed score squared components
### hessian_list: list of all precomputed hessian components 
### indices: list of indices to compute for 
### q: # of phenotypes
build_variance_matrix <- function(ss_list, hessian_list, indices, q) {
  ### Compute each v_ij entry we need to build out the matrix
  indices$v_ij <- 
    purrr::map_dbl(indices$row_id, ~{
      i <- indices$i[.x]
      j <- indices$j[.x]
      
      compute_variance_entry(score_squared = ss_list[[.x]],
                             hessian_1 = hessian_list[[i]],
                             hessian_2 = hessian_list[[j]])
    })
  
  ### Build Matrix
  V <- matrix(0, nrow = q, ncol = q)
  V[lower.tri(V, diag = T)] <- indices$v_ij
  V <- V + t(V)
  diag(V) <- diag(V)/2
  
  return(V)
}
