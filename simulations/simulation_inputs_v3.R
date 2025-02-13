library(readr)
library(purrr)

build_beta <- function(beta_seq, n, start = 1, by = 2) {
  beta_mat <- NULL
  for(i in 1:length(beta_seq)) {
    m <- matrix(unlist(map(seq(start,n,by), ~rep(c(beta_seq[i], 0), c(.x, n-.x)))),
                nrow = ceiling(n/by), byrow = T)
    beta_mat <- rbind(beta_mat, m)
  }
  return(beta_mat)
}



### Simulation 1a: Binary Phenotypes and Continuous Phenotypes
### Same Dir Sims
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 4,

               'num_pca' = rep(4, 5),
               'na_freq' = 0.10,
               'maf' = 0.20,
               'prevalence' = 0.30,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = F,
               'missingness' = 'mcar')

grain <- 40
beta <- build_beta(seq(0.01, 0.35, length.out = grain), n = 4, by = 1, start = 1)
ix_1 <- ( apply(beta, 1, function(x) sum(x != 0)) == 1 )
ix_2 <- ( apply(beta, 1, function(x) sum(x != 0)) == 2 )
ix_3 <- ( apply(beta, 1, function(x) sum(x != 0)) == 3 )
beta_opp <- beta
beta_opp[, c(2,4)]  <- beta_opp[, c(2,4)] * -1
beta_opp[ix_1, 1] <- -1 * beta_opp[ix_1, 1]
beta_con <- beta[!ix_3,]
beta_bin <- beta[!ix_2,]

Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- 0.4
Sigma[5:8, 5:8] <- 0.4
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_1_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_1_inputs, 'inputs/v3/simulation_run_v3_same_sign_1.rds')

### No Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_2_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_2_inputs, 'inputs/v3/simulation_run_v3_same_sign_2.rds')

### Negative Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- -0.3
Sigma[5:8, 5:8] <- -0.3
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_3_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_3_inputs, 'inputs/v3/simulation_run_v3_same_sign_3.rds')


### Simulation 1b: Binary Phenotypes and Continuous Phenotypes
### Opp Dir Sims
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 4,

               'num_pca' = rep(4, 5),
               'na_freq' = 0.10,
               'maf' = 0.20,
               'prevalence' = 0.30,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = F,
               'missingness' = 'mcar')

grain <- 40
beta <- build_beta(seq(0.01, 0.35, length.out = grain), n = 4, by = 1, start = 1)
ix_1 <- ( apply(beta, 1, function(x) sum(x != 0)) == 1 )
ix_2 <- ( apply(beta, 1, function(x) sum(x != 0)) == 2 )
ix_3 <- ( apply(beta, 1, function(x) sum(x != 0)) == 3 )
beta_opp <- beta
beta_opp[, c(2,4)]  <- beta_opp[, c(2,4)] * -1
beta_opp[ix_1, 1] <- -1 * beta_opp[ix_1, 1]
beta_con <- beta_opp[!ix_3,]
beta_bin <- beta[!ix_2,]

Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- 0.4
Sigma[5:8, 5:8] <- 0 ### No correlation for binary
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_opp_sign_1_inputs <- list('params' = params,
                              'beta_bin' = beta_bin,
                              'beta_con' = beta_con)

write_rds(sim_opp_sign_1_inputs, 'inputs/v3/simulation_run_v3_opp_sign_1.rds')

### No Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_opp_sign_2_inputs <- list('params' = params,
                              'beta_bin' = beta_bin,
                              'beta_con' = beta_con)

write_rds(sim_opp_sign_2_inputs, 'inputs/v3/simulation_run_v3_opp_sign_2.rds')

### Negative Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- -0.3
Sigma[5:8, 5:8] <- 0 ### No correlation for binary
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_opp_sign_3_inputs <- list('params' = params,
                              'beta_bin' = beta_bin,
                              'beta_con' = beta_con)

write_rds(sim_opp_sign_3_inputs, 'inputs/v3/simulation_run_v3_opp_sign_3.rds')



### Simulation 2: Binary Phenotypes Only, Common Variant
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 8,

               'num_pca' = rep(4, 5),
               'na_freq' = 0.10,
               'maf' = 0.2,
               'prevalence' = 0.3,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = F,
               'missingness' = 'mcar')

grain <- 30
beta <- build_beta(seq(0.01, 0.6, length.out = grain), n = 8, by = 3, start = 2)
beta_neg <- beta
beta_neg[, seq(2, 8, 3)] <- -beta_neg[, seq(2, 8, 3)]
beta_all_neg <- -beta
zero_mat <- build_beta(rep(0, grain), n = 8, by = 3, start = 2)
beta_con <- rbind(zero_mat, zero_mat, zero_mat)
beta_bin <- rbind(beta, beta_all_neg, beta_neg)

Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- 0.4
Sigma[5:8, 5:8] <- 0.7
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_2_inputs <- list('params' = params,
                     'beta_bin' = beta_bin,
                     'beta_con' = beta_con)

write_rds(sim_2_inputs, 'inputs/v3/simulation_run_v3_2.rds')

### Simulation 3: Binary Phenotypes Only, Rare Variant
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 8,

               'num_pca' = rep(4, 5),
               'na_freq' = 0.10,
               'maf' = 0.05,
               'prevalence' = 0.10,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = F,
               'missingness' = 'mcar')

grain <- 30
beta <- build_beta(seq(0.01, 0.6, length.out = grain), n = 8, by = 3, start = 2)
beta_neg <- beta
beta_neg[, seq(2, 8, 3)] <- -beta_neg[, seq(2, 8, 3)]
beta_all_neg <- -beta
zero_mat <- build_beta(rep(0, grain), n = 8, by = 3, start = 2)
beta_con <- rbind(zero_mat, zero_mat, zero_mat)
beta_bin <- rbind(beta, beta_all_neg, beta_neg)

Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- 0.4
Sigma[5:8, 5:8] <- 0.7
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_3_inputs <- list('params' = params,
                     'beta_bin' = beta_bin,
                     'beta_con' = beta_con)

write_rds(sim_3_inputs, 'inputs/v3/simulation_run_v3_3.rds')

### Simulation 4: Binary Phenotypes Only, Healthy Controls
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 8,

               'num_pca' = rep(4, 5),
               'na_freq' = 0,
               'maf' = 0.2,
               'prevalence' = 0.30,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = T,
               'missingness' = 'mcar')

grain <- 30
beta <- build_beta(seq(0.01, 0.6, length.out = grain), n = 8, by = 3, start = 2)
beta_neg <- beta
beta_neg[, seq(2, 8, 3)] <- -beta_neg[, seq(2, 8, 3)]
beta_all_neg <- -beta
zero_mat <- build_beta(rep(0, grain), n = 8, by = 2, start = 2)
beta_con <- rbind(zero_mat, zero_mat, zero_mat)
beta_bin <- rbind(beta, beta_all_neg, beta_neg)

Sigma <- matrix(0.7, nrow = 8, ncol = 8)
Sigma[1:4, 5:8] <- 0.4
Sigma[5:8, 1:4] <- 0.4
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_4_inputs <- list('params' = params,
                     'beta_bin' = beta_bin,
                     'beta_con' = beta_con)

write_rds(sim_4_inputs, 'inputs/v3/simulation_run_v3_4.rds')

### Simulation 5: Re-run Sim 1 but MAR rather than MCAR
### Same Dir Sims
params <- list('num_sites' = 5,
               'n_k' = rep(1000, 5),

               'q' = 8,
               'q_k' = c(8,7,6,5,8),
               'q_bin' = 4,

               'num_pca' = rep(4, 5),
               'na_freq' = 0.10,
               'maf' = 0.20,
               'prevalence' = 0.30,

               'rm_gender' = F,
               'flip_y' = F,
               'intercept_only' = F,
               'healthy_controls' = F,
               'missingness' = 'mar')

grain <- 40
beta <- build_beta(seq(0.01, 0.35, length.out = grain), n = 4, by = 1, start = 1)
ix_1 <- ( apply(beta, 1, function(x) sum(x != 0)) == 1 )
ix_2 <- ( apply(beta, 1, function(x) sum(x != 0)) == 2 )
ix_3 <- ( apply(beta, 1, function(x) sum(x != 0)) == 3 )
beta_opp <- beta
beta_opp[, c(2,4)]  <- beta_opp[, c(2,4)] * -1
beta_opp[ix_1, 1] <- -1 * beta_opp[ix_1, 1]
beta_con <- beta[!ix_3,]
beta_bin <- beta[!ix_2,]

Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- 0.4
Sigma[5:8, 5:8] <- 0.4
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_1_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_1_inputs, 'inputs/v3/simulation_run_v3_mar_1.rds')

### No Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_2_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_2_inputs, 'inputs/v3/simulation_run_v3_mar_2.rds')

### Negative Correlation
Sigma <- matrix(0, nrow = 8, ncol = 8)
Sigma[1:4, 1:4] <- -0.3
Sigma[5:8, 5:8] <- -0.3
diag(Sigma) <- 1
params$Sigma <- Sigma

sim_same_sign_3_inputs <- list('params' = params,
                               'beta_bin' = beta_bin,
                               'beta_con' = beta_con)

write_rds(sim_same_sign_3_inputs, 'inputs/v3/simulation_run_v3_mar_3.rds')
