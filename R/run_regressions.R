### Functions to Run Reduced Models

### Run a Single Regression for one Phenotype
###
### x: vector of SNPs in [0,1,2]
### y: vector of phenotypes
### Z: matrix of site-specific, pheno-type specific covariates
### delta: vector indicating if phenotype is present (1) or missing (0)
### type: specifies the type of phenotype (e.g. continuous, binary, count)

single_regression <- function(x, y, Z, delta, type) {
  ncz <- ncol(Z)
  ### Fit model according to the type
  if(type == 'continuous') { ### Intercept Only
    if(ncz == 0) {
      model <- lm(y[delta == 1] ~ 1)
    } else {
      model <- lm(y[delta == 1] ~ Z[delta == 1,])
    }
  } else if(type == 'binary') {
    if(ncz == 0) { ### Intercept Only
      model <- glm(y[delta == 1] ~ 1, family = binomial(link = 'logit'))
    } else {
      model <- glm(y[delta == 1] ~ Z[delta == 1,], family = binomial(link = 'logit'))
    }
  } else if(type == 'count') {
    if(ncz == 0) { ### Intercept Only
      model <- glm(y[delta == 1] ~ 1, family = poisson(link = 'log'))
    } else {
      model <- glm(y[delta == 1] ~ Z[delta == 1,], family = poisson(link = 'log'))
    }
  }

  ### Extract Coefficients
  gamma <- model$coefficients

  return(gamma)
}

### Run all regressions at once
###
### x: vector of SNPs in [0,1,2]
### Y: Matrix of phenotypes (one per column)
### Z_list: list of matrices, each corresponding to a set of site-specific, pheno-type specific covariates (
### Z_index: index indicating which set of covariates can be shared by multiple phenotypes
### Delta: vector indicating if phenotype is present (1) or missing (0)
### types: vector specifies the type of phenotype (e.g. continuous, binary, count)

run_regressions <- function(x, Y, Z_list, Z_index, Delta, types) {
  gamma_list <-
    purrr::map(1:ncol(Y), ~single_regression(x = x,
                                             y = Y[,.x],
                                             Z = Z_list[[Z_index[.x]]],
                                             delta = Delta[,.x],
                                             type = types[.x]))

  return(gamma_list)
}
