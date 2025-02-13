source('helpers.R')
library(ASSET)
library(harmonicmeanp)

pheWAS <- function(df_phenotypes, df_covariates, col_types, site = 'mega') {
  if(!is.data.frame(df_phenotypes)) {
    df_phenotypes <- as_tibble(df_phenotypes)
  }
  if(!is.data.frame(df_covariates)) {
    df_covariates <- as_tibble(df_covariates)
  }

  df_phewas <- map_dfr(1:ncol(df_phenotypes), ~single_regression(y = df_phenotypes[[.x]],
                                                                 covariates = df_covariates,
                                                                 type = col_types[.x],
                                                                 site = site))
  df_phewas$site <- site
  df_phewas$phenotype <- colnames(df_phenotypes)

  return(df_phewas)
}


run_pheWAS <- function(snps, phenotypes, covariates, n_k, run_asset = F, return_z = F) {
  ### PheWAS mega
  df_phenotypes <-
    map_dfr(phenotypes, as_tibble)

  df_covariates <-
    map_dfr(covariates, as_tibble) %>%
    mutate('snp' = unlist(snps)) %>%
    mutate('site_id' = as.character(rep(1:length(n_k), n_k)))

  col_types <- unlist(infer_types(df_phenotypes))

  pheWAS_mega <- pheWAS(df_phenotypes, df_covariates, col_types)
  pheWAS_mega$z <- pheWAS_mega$beta/sqrt(pheWAS_mega$var)

  ### Min P for Phewas Mega
  min_p_mega <- min(pheWAS_mega$p)
  q <- length(pheWAS_mega$p)

  ### ACAT of P-Values
  acat_p_mega <- acat(pheWAS_mega$p)

  ### GHC/GBJ
  gbj_mega <- GBJ::GBJ(pheWAS_mega$z, diag(1, length(pheWAS_mega$z)))$GBJ_pvalue
  ghc_mega <- GBJ::GHC(pheWAS_mega$z, diag(1, length(pheWAS_mega$z)))$GHC_pvalue

  ### ACAT Omnibus Test
  acat_p_pheWAS_mega <- acat(c(min(1, min_p_mega * q), gbj_mega, ghc_mega, acat_p_mega))

  ### PheWAS meta
  col_types <- infer_types(phenotypes)
  covariates <- map2(covariates, snps, ~cbind(.x, 'snp' = .y))

  pheWAS_meta <-
    map_dfr(1:length(phenotypes), ~pheWAS(phenotypes[[.x]],
                                          covariates[[.x]],
                                          col_types[[.x]],
                                          site = .x))

  pheWAS_meta <-
    pheWAS_meta %>%
    mutate('weight' = 1/var) %>%
    group_by(phenotype) %>%
    summarise('beta_meta' =  weighted.mean(beta, weight),
              'var_meta' = 1/sum(weight)) %>%
    mutate('z' = beta_meta/sqrt(var_meta)) %>%
    mutate('p' = 2 * pnorm(-abs(z)))

  ### Min P
  min_p_meta <- min(pheWAS_meta$p)

  ### ACAT of P-Values
  acat_p_meta <- acat(pheWAS_meta$p)

  ### GHC/GBJ
  gbj_meta <- GBJ::GBJ(pheWAS_meta$z, diag(1, length(pheWAS_meta$z)))$GBJ_pvalue
  ghc_meta <- GBJ::GHC(pheWAS_meta$z, diag(1, length(pheWAS_meta$z)))$GHC_pvalue

  ### ACAT Omnibus test
  acat_p_pheWAS_meta <- acat(c(min(1, min_p_meta * q), gbj_meta, ghc_meta, acat_p_meta))

  ### ASSET (Subset Based Meta Analysis Test)
  if(run_asset) {
    asset_p_meta <- asset(phenotypes, pheWAS_meta, q)
  } else {
    asset_p_meta <- NA
  }

  ### Harmonic Min P
  hmp <- hmp.stat(pheWAS_meta$p)


  if(!return_z) {
    return(list('p_pheWAS_mega' = min_p_mega,
                'p_pheWAS_mega_acat' = acat_p_pheWAS_mega,
                'p_pheWAS_mega_acat_p' = acat_p_mega,

                'p_pheWAS_meta' = min_p_meta,
                'p_pheWAS_meta_acat' = acat_p_pheWAS_meta,
                'p_pheWAS_meta_acat_p' = acat_p_meta,

                'p_asset_meta' = asset_p_meta,
                'p_pheWAS_hmp' = hmp))
  } else {
    df <-
      pheWAS_mega %>%
      select(phenotype, 'mega' = z) %>%
      inner_join(select(pheWAS_meta, phenotype, 'meta' = z), by = 'phenotype')
    return(df)
  }
}
