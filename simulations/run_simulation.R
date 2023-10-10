library(furrr)
library(mixWAS)
library(dplyr)
library(ASSET)
source('pheWAS.R')
source('generate_data.R')

run_simulation <- function(params, beta_bin, beta_con, n_sims, alpha = 0.05, verbose = F, future_seeds) {
  
  cat(ifelse(verbose, 'Generating Datasets\n', ''))
  ### Generate Datasets 
  datasets <- future_map(1:n_sims, ~generate_data(num_sites = params$num_sites, 
                                                  n_k = params$n_k, 
                                                  q = params$q,
                                                  q_k = params$q_k,
                                                  q_bin = params$q_bin, 
                                                  maf = params$maf,
                                                  prevalence = params$prevalence, 
                                                  na_freq = params$na_freq,
                                                  num_pca = params$num_pca,
                                                  
                                                  rm_gender = params$rm_gender,
                                                  intercept_only = params$intercept_only,
                                                  flip_y = params$flip_y,
                                                  
                                                  beta_bin = beta_bin,
                                                  beta_con = beta_con,
                                                  Sigma = params$Sigma,
                                                  
                                                  healthy_controls = params$healthy_controls,
                                                  
                                                  sim_id = .x), 
                         .options = furrr_options(seed = future_seeds[1]))
  
  
  ### Empirical Prevalence
  df_ep <- 
    map_dfr(datasets, ~pluck(.x, 'df_ep')) %>% 
    summarise(across(everything(), ~mean(.x, na.rm = T)))
  
  cat(ifelse(verbose, 'Running mixWAS Component Return Mode\n', ''))
  
  ### MixWAS Components (Easier to Parallelize prior to p-values)
  components <- future_map(datasets, ~mixWAS(snps = .x$snps,
                                             phenotypes = .x$phenotypes,
                                             covariates = .x$covariates, 
                                             return_p = F))
  
  cat(ifelse(verbose, 'Computing p-values\n', ''))
  ### P-value computation
  p_values <- map_dfr(components, ~run_hypothesis_tests(.x$score, .x$V_inv, .x$z, .x$q))
  
  ### PheWAS
  cat(ifelse(verbose, 'Running PheWAS\n', ''))
  phewas <- 
    future_map_dfr(datasets, ~run_pheWAS(snps = .x$snps,
                                         phenotypes = .x$phenotypes,
                                         covariates = .x$covariates,
                                         run_asset = params$q_bin == params$q,
                                         n_k = params$n_k), 
                   .options = furrr_options(seed = future_seeds[2])) 
  p_values <- bind_cols(p_values, phewas)
  
  ### Oracle Test
  cat(ifelse(verbose, 'Running Oracle test\n', ''))
  if(params$q == params$q_bin) {
    ix_non_null <- beta_bin != 0
  } else if(params$q > params$q_bin & params$q_bin > 0) {
    ix_non_null <- c(beta_con != 0, beta_bin != 0)
  } else {
    ix_non_null <- beta_bin != 0
  }
  
  p_oracle <- future_map_dbl(components, ~oracle_test(score = .x$score,
                                                      V_inv = .x$V_inv,
                                                      ix_non_null = ix_non_null))
  
  p_oracle_uncorrelated <- future_map_dbl(components, ~oracle_test_uncorrelated(z_score = .x$z,
                                                                                ix_non_null = ix_non_null))
  
  
  p_values$p_oracle <- p_oracle
  p_values$p_oracle_uncorrelated <- p_oracle_uncorrelated ### Decorrelate Oracle P-value
  
  ### ASSET Subset based test (score based version)
  # cat(ifelse(verbose, 'Running ASSET test\n', ''))
  # p_asset_score <- future_map2_dbl(components, datasets, asset_score) ### Score based version
  # p_values$p_asset_score <- p_asset_score
  
  
  ### Compoling Info About Simulation
  cat(ifelse(verbose, 'Saving Sim Info\n', ''))  
  df_info <- 
    tibble('alpha' = alpha, 
           'n_sims' = n_sims, 
           
           'num_sites' = params$num_sites,
           'total_subjects' = sum(params$n_k),
           'max_site' = max(params$n_k),
           
           'q' = params$q,
           'q_bin' = params$q_bin,
           'q_con' = params$q - params$q_bin,
           
           'intercept_only' = params$intercept_only,
           'flip_y' = params$flip_y,
           'rm_gender' = params$rm_gender,
           
           'na_freq' = params$na_freq,
           'maf' = params$maf,
           'prevalence' = params$prevalence,
           'num_pca' = params$num_pca[1],
           
           'n_true_con' = sum(beta_con != 0),
           'max_beta_con' = max(abs(beta_con)),
           'direction_con' = case_when(any(sign(beta_con) > 0) & any(sign(beta_con) < 0) ~ 'Opposite', 
                                       all(sign(beta_con) >= 0) ~ 'Same -- Positive',
                                       T ~ 'Same -- Negative'),
           
           
           
           'n_true_bin' = sum(beta_bin != 0),
           'max_beta_bin' = max(abs(beta_bin)),
           'direction_bin' = case_when(any(sign(beta_bin) > 0) & any(sign(beta_bin) < 0) ~ 'Opposite', 
                                       all(sign(beta_bin) >= 0) ~ 'Same -- Positive',
                                       T ~ 'Same -- Negative'),
           
           'beta_con' = paste(beta_con, collapse = '/'),
           'beta_bin' = paste(beta_bin, collapse = '/'),
           
           'n_k' = paste(params$n_k, collapse = '/'),
           'q_k' = paste(params$q_k, collapse = '/'))
  
  ### Compute Power
  cat(ifelse(verbose, 'Computing Power\n', ''))
  df_power <- 
    p_values %>%
    summarise(across(.cols = -c('p_pheWAS_mega', 'p_pheWAS_meta'), .fns = ~mean(.x <= alpha))) %>% 
    ### Bonferonni Correction for PheWAS min_p methods
    bind_cols(summarise(p_values, across(.cols = c('p_pheWAS_mega', 'p_pheWAS_meta'), .fns = ~mean(.x <= alpha/params$q)))) %>% 
    bind_cols(df_info) %>% 
    bind_cols(df_ep)
  
  
  return(df_power)
}
