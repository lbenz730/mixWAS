### Example to show how mixWAS works
### Generate a dataset for illustration  

source('mixWAS.R')
library(magrittr)

# future::plan(future::multicore(workers = 15))
# options(future.fork.enable = T)
set.seed(123)

### Parameters
# num_sites <- 2
# n_k <- c(100, 100)
# maf <- 0.2
# q <- 5
# q_k <- c(5, 3)
# q_bin <- 3
# q_con <- q - q_bin
# prevelance <- runif(q_bin, min = 0.1, 0.3)
# na_freq <- 0.01
# n_covariates <- 5

# num_sites <- 2
# n_k <- c(10000, 10000)
# maf <- 0.2
# q <- 5
# q_k <- c(5, 3)
# q_bin <- 3
# q_con <- q - q_bin
# prevelance <- runif(q_bin, min = 0.1, 0.3)
# na_freq <- 0.01
# n_covariates <- 5

# num_sites <- 2
# n_k <- c(10000, 10000)
# 
# maf <- 0.2
# q <- 20
# q_k <- c(20, 15)
# q_bin <- 5
# q_con <- q - q_bin
# prevelance <- runif(q_bin, min = 0.1, 0.3)
# na_freq <- 0.01
# n_covariates <- 10

# num_sites <- 5
# n_k <- c(5000, 250, 1700, 3000, 10000)
# maf <- 0.2
# q <- 20
# q_k <- c(20, 10, 17, 9, 12)
# q_bin <- 10
# q_con <- q - q_bin
# prevelance <- runif(q_bin, min = 0.1, 0.3)
# na_freq <- 0.01
# n_covariates <- 5


### Generate SNPs
snps <- purrr::map(n_k, ~{rbinom(.x, 1, 1-maf) + rbinom(.x, 1, 1-maf)})

### Phenotypes 
binary <- 
  purrr::map(n_k, function(n){
    purrr::map2_dfc(1:q_bin, prevelance, ~rbinom(n, 1, .y)) %>% 
      set_names(paste0('Y', q_con + 1:q_bin))
  })

continuous <- 
  purrr::map(n_k, function(n) {
    purrr::map_dfc(1:q_con, ~rnorm(n)) %>% 
      set_names(paste0('Y', 1:q_con))
  })

phenotypes <- purrr::map(1:num_sites, ~{
  df <- dplyr::bind_cols(continuous[[.x]], binary[[.x]])
  ### Randomly Have a few of the phenotypes be missing 
  na_ind <- matrix(runif(nrow(df) * ncol(df)), ncol = ncol(df)) < na_freq
  df[na_ind] <- NA
  # Randomly choose q_k phenotypes in the site
  # df <- df[,sample(1:q, q_k[.x])]
  as.matrix(df)
})

### Covariates
covariates <- 
  purrr::map(n_k, function(n) {
    purrr::map_dfc(1:n_covariates, ~rnorm(n)) %>% 
      set_names(paste0('Z', 1:n_covariates)) %>% 
      as.matrix()
  })

# ### Make list of a few covariates for each phenotype
# covariates <- 
#   purrr::map2(covariates, q_k, function(M, qk) { 
#     purrr::map(1:qk, ~M[, sample(1:ncol(M), sample(2:ncol(M), 1))])
#   })

out_new <- 
  mixWAS_single_site(snps = snps[[1]],
                     phenotypes = phenotypes[[1]],
                     covariates = covariates[[1]])

# mixWAS_single_site(snps = snps[[2]],
#                    phenotypes = phenotypes[[2]],
#                    covariates = covariates[[2]])
# 
# t0 <- Sys.time()
# mixWAS(snps, phenotypes, covariates)                    
# t1 <- Sys.time()                  
# t1-t0


X_site1 <- snps[[1]]
Y_site1 <- phenotypes[[1]]
Delta_site1 <- 1 - is.na(Y_site1) 
Y_site1[is.na(Y_site1)] <- 999
Zlist_site1 <- purrr::map(1:q, ~covariates[[1]])

out_old <- run.each.site(X_site1,Y_site1,Zlist_site1,Delta_site1,q1=q_con)

out_new$score - out_old[[1]]
all(abs(out_new$V - out_old[[2]]) < 1e-12)
