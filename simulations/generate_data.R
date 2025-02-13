library(dplyr)
library(tidyr)
library(purrr)
library(bindata)
source('helpers.R')

expit <- function(p) {
  exp(p)/(exp(p) + 1)
}

logit <- function(p) {
  log(p/(1-p))
}

### Function to generate data to use in mixWAS power simulations
###
### num_sites: number of sites
### n_k: vector of length num_sites
### q: # of phenotypes
### q_k: Vector of phenotypes to choose for each site (chose q_k randomly for site k)
### q_bin: # of binomial covariates (q_bin + q_con = q)
### maf: minor allele frequency
### na_freq: % of covariates to randomly set to missing
### gamma_bin: Binary Patient Covariates (list, for each site/phenotype)
### gamma_con: Continuous Patient Covariates (list, for each site/phenotype)
### num_pca: Number of Principle components
### beta_bin: beta vector for binary covariates
###       and the remaining q_bin are for binary
###
### New Diagnostic arguments for debugging rare variant cases (2022-03-11)
### flip_y: Y = 1 - Y
### intercept_only: No Z
### rm_gender: remove gender
###
### Sigma: Optional correlation matrix for rmvbin()/mvrnorm
###
### New Healthy Controls Argument (2022-07-07)
### healthy_controls: Flag to force us to take controls to be only patients
###                   that are controls (0) for all phenotypes. If TRUE, we
###                   we will NA out any 0 in the phenotypes for patients with
###                   case status (1) for at least one of the phenotypes
###
### Missingness (default = 'mcar', alternative option of 'mar')



generate_data <- function(num_sites, n_k,
                          q, q_k, q_bin,
                          maf, prevalence, na_freq,
                          num_pca, beta_bin, beta_con,
                          sim_id,
                          flip_y = F,
                          intercept_only = F,
                          rm_gender = F,
                          Sigma = NULL,
                          healthy_controls = F,
                          missingness = 'mcar') {

  ### Number of continuous covariates
  q_con <- q - q_bin

  if(q_con == 0) {
    n_pos <- sum(beta_bin > 0)
    n_neg <- sum(beta_bin < 0)

    if(n_neg != 0 & n_pos != 0 & n_neg != n_pos) {
      beta_bin[n_neg + n_pos] <- sample(c(-1, 1), 1) * beta_bin[n_neg + n_pos]
    }
  }

  ### Some Random Gammas
  gammas_bin <-
    map(num_pca, function(n_var){ map(1:q_bin, ~{
      gamma <-
        c(runif(n_var, -0.5, 0.5), ### PCs
          runif(1, -0.05, 0.05)) ### Age

      if(!rm_gender) {
        gamma <- c(gamma, runif(1, -0.1, 0.1)) ### gender
      }

      gamma

    }) })

  gammas_con <-
    map(num_pca, function(n_var){ map(1:q_con, ~{
      gamma <-
        c(runif(n_var, -0.5, 0.5), ### PCs
          runif(1, -0.05, 0.05)) ### Age

      if(!rm_gender) {
        gamma <- c(gamma, runif(1, -0.1, 0.1)) ### gender
      }

      gamma

    }) })


  ### Control Flags for diagnostic testing
  rm_vars <- c()
  if(intercept_only) {
    gammas_con <- list()
    gammas_bin <- list()
    rm_vars <- c(paste0('PC', 1:max(num_pca)), 'age_centered', 'gender')
  } else if(rm_gender) {
    rm_vars <- 'gender'
  }


  ### Generate SNPs from binomial(2, maf)
  ### Center by subracting E(SNP)
  snps <- map(n_k, ~{rbinom(.x, 2, maf) - 2 * maf})


  ### Covariates
  covariates <-
    map2(n_k, num_pca, function(n, nc) {
      map_dfc(1:nc, ~{
        df <-
          tibble('z' = rnorm(n)) %>%
          set_names(paste0('PC', .x))
      }) %>%
        mutate('age_centered' = rnorm(n, mean = 0, sd = 15), ### Centered age
               'gender' =  rbinom(n, 1, 0.5)) %>%
        select(-any_of(rm_vars)) %>% ### remove gender or all variables
        as.matrix()
    })

  ### Phenotypes
  phenotypes <- list()
  emp_prev <- c()
  no_cases <- c()
  for(site in 1:num_sites) {
    m <- matrix(NA, nrow = n_k[site], ncol = q)

    if(is.null(Sigma)) {
      ### Continuous Phenotypes
      if(q_con > 0) {
        for(phenotype in 1:q_con) {
          if(!intercept_only) {
            m[,phenotype] <- covariates[[site]] %*% gammas_con[[site]][[phenotype]] + snps[[site]] * beta_con[phenotype] + rnorm(n_k[site], mean = 0, sd = 2.3)
          } else { ### Intercept Only
            m[,phenotype] <- snps[[site]] * beta_con[phenotype] + rnorm(n_k[site], mean = 0, sd = 2.3)
          }
        }
      }

      ### Binary Phenotypes
      if(q_bin > 0) {
        for(phenotype in 1:q_bin) {
          if(!intercept_only) {
            log_odds <- covariates[[site]] %*% gammas_bin[[site]][[phenotype]] + snps[[site]] * beta_bin[phenotype] + logit(prevalence)

          } else { ### Intercept Only
            log_odds <- snps[[site]] * beta_bin[phenotype] + logit(prevalence)
          }
          p <- expit(log_odds)
          m[,phenotype + q_con] <- rbinom(n = length(p), size = 1, prob = p)
        }
      }
    } else { ### Correlated Phenotypes
      if(q_con > 0) { ### Continuous
        for(phenotype in 1:q_con) {
          if(!intercept_only) {
            m[,phenotype] <- covariates[[site]] %*% gammas_con[[site]][[phenotype]] + snps[[site]] * beta_con[phenotype]
          } else { ### Intercept Only
            m[,phenotype] <- snps[[site]] * beta_con[phenotype]
          }
        }

        m[,1:q_con] <- m[,1:q_con] + MASS::mvrnorm(n = n_k[site], mu = rep(0, q_con), Sigma = Sigma[1:q_con, 1:q_con] * 2.3^2)

      }
      if(q_bin > 0) {
        ### Correlated Phenotypes (binary)
        prob <- matrix(NA, nrow = n_k[site], ncol = q_bin)

        for(phenotype in 1:q_bin) {
          if(!intercept_only) {
            log_odds <- covariates[[site]] %*% gammas_bin[[site]][[phenotype]] + snps[[site]] * beta_bin[phenotype] + logit(prevalence)

          } else { ### Intercept Only
            log_odds <- snps[[site]] * beta_bin[phenotype] + logit(prevalence)
          }
          prob[, phenotype] <- expit(log_odds)
        }


        ### Sample Phenotypes
        m[,(q_con + 1):q] <- t(apply(prob, 1, function(mu) { rmvbin(1, margprob = mu, sigma = Sigma[(q_con+1):q, (q_con+1):q] ) } ))
      }
    }

    colnames(m) <- paste0('Y', 1:q)


    ### Measure the Empirical Prevalence
    emp_prev[site] <- mean(m[,q_con + 1:q_bin], na.rm = T)
    no_cases[site] <- mean(apply(as.matrix(m[,q_con + 1:q_bin]), 2, function(x){ all(x == 0) | all(x == 1) }))

    ### Randomly choose q_k phenotypes in the site
    if(q_k[site] > 1) {
      m <- m[,sort(sample(1:q, q_k[site]))]
    }

    ### Randomly Have a few of the phenotypes be missing
    if(missingness == 'mcar') {
      na_ind <- matrix(runif(nrow(m) * ncol(m)), ncol = ncol(m)) < na_freq
      m[na_ind] <- NA
    } else if(missingness == 'mar') {
      p_miss <- expit( logit(na_freq * 0.8) + 0.02 * covariates[[site]][,'age_centered'] )
      na_ind <- matrix(runif(nrow(m) * ncol(m)), ncol = ncol(m)) < matrix(rep(p_miss, ncol(m)), ncol = ncol(m), byrow = F)
      m[na_ind] <- NA
    }

    ### Diagnostic to flip Y
    if(flip_y) {
      m <- 1 - m
    }

    phenotypes[[site]] <- m

  }

  ### Apply Healthy Controls
  if(healthy_controls) {
    phenotypes <- get_healthy_controls(phenotypes)
  }

  ### Empirical prevalence
  df_ep <- tibble('empirical_prev' = weighted.mean(emp_prev, n_k),
                  'no_case_rate' = weighted.mean(no_cases, q_k))

  return(list('snps' = snps,
              'covariates' = covariates,
              'phenotypes' = phenotypes,
              'sim_id' = sim_id,
              'df_ep' = df_ep))

}
