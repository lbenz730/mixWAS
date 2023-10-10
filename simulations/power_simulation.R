library(readr)
library(furrr)
library(dplyr)
library(purrr)

set.seed(123)
n_cores <- 32
plan(future::multicore(workers = n_cores))
options(future.fork.enable = T)
options(future.globals.maxSize = 8 * 1024^3)

source('generate_data.R')
source('run_simulation.R')

args <- commandArgs(trailingOnly = T)
sim_id <- as.numeric(args[1])
sim_type <- args[2]
if(is.na(sim_type)) {
  simulation_run <- paste0('v3_', sim_id)
} else {
  simulation_run <- paste0('v3_', sim_type, '_', sim_id)
}
inputs <- read_rds(paste0('inputs/v3/simulation_run_', simulation_run, '.rds'))
params <- inputs$params
beta_bin <- inputs$beta_bin
beta_con <- inputs$beta_con
n <- nrow(beta_bin)

#### Check if Job was Requeued
if(dir.exists(paste0('outputs/v3/simulation_run_', simulation_run))) {
  files <- dir(paste0('outputs/v3/simulation_run_', simulation_run))
  start <- max(as.numeric(gsub('[^0-9]', '', files[!grepl('combined', files)])), na.rm = T)
} else {
  start <- 0 
}

for(i in 1:n) {
  cat('Running Simulation', i, 'of', n, '\n') 
  ### Seed to use in parallel
  seeds <- sample(1:1e7, 2)
  
  if(i > start) {
    
    df <- 
      run_simulation(params, 
                     beta_bin = beta_bin[i,],
                     beta_con = beta_con[i,],
                     n_sims = 2000, 
                     verbose = T, 
                     future_seeds = seeds)
    
    if(!dir.exists(paste0('outputs/v3/simulation_run_', simulation_run))) {
      dir.create(paste0('outputs/v3/simulation_run_', simulation_run))
    }
    write_csv(df, paste0('outputs/v3/simulation_run_', simulation_run, '/', i, '.csv'))
    
    
    ## Clean up memory
    gc()
    
  }
}

### Save Combined File
files <- dir(paste0('outputs/v3/simulation_run_', simulation_run), full.names = T)
files <- files[!grepl('combined', files)]
df_combined <-
  map_dfr(files, read_csv)
write_csv(distinct(df_combined), paste0('outputs/v3/simulation_run_', simulation_run, '/combined.csv'))
