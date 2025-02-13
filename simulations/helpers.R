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

single_regression <- function(y, covariates, type, site = 'mega') {
  df <- covariates
  df$y <- y

  df <- filter(df, !is.na(y))
  if(!'site_id' %in% names(df)) {
    df$site_id <- '1'
  }

  ### Fit Model
  if(site == 'mega' & n_distinct(df$site_id) > 1 & 'gender' %in% names(df)) {
    if(type == 'continuous') {
      model <- lm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + gender:site_id + snp, data = df)
    } else if(type == 'binary') {
      model <- glm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + gender:site_id + snp, data = df, family = binomial(link = 'logit'))
    } else if(type == 'count') {
      model <- glm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + gender:site_id + snp, data = df, family = poisson(link = 'log'))
    }
  } else if(site == 'mega' & n_distinct(df$site_id) > 1 & 'PC1' %in% names(df)) {
    if(type == 'continuous') {
      model <- lm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + snp, data = df)
    } else if(type == 'binary') {
      model <- glm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + snp, data = df, family = binomial(link = 'logit'))
    } else if(type == 'count') {
      model <- glm(y ~ site_id + PC1:site_id + PC2:site_id + PC3:site_id + PC4:site_id + age_centered:site_id + snp, data = df, family = poisson(link = 'log'))
    }
  } else if(site == 'mega' & n_distinct(df$site_id) > 1) {
    if(type == 'continuous') {
      model <- lm(y ~ site_id + snp, data = df)
    } else if(type == 'binary') {
      model <- glm(y ~ site_id + snp, data = df, family = binomial(link = 'logit'))
    } else if(type == 'count') {
      model <- glm(y ~ site_id + snp, data = df, family = poisson(link = 'log'))
    }
  } else{
    df <- select(df, -any_of(c('site_id')))
    if(type == 'continuous') {
      model <- lm(y ~ ., data = df)
    } else if(type == 'binary') {
      model <- glm(y ~., data = df, family = binomial(link = 'logit'))
    } else if(type == 'count') {
      model <- glm(y ~ ., data = df, family = poisson(link = 'log'))
    }
  }


  ### Extract Coefficients
  model_summary <- summary(model)$coefficients

  beta_snp <- model_summary['snp', 'Estimate']
  var_snp <- model_summary['snp', 'Std. Error']^2
  if(type == 'continuous') {
    p_snp <- model_summary['snp', 'Pr(>|t|)']
  } else if(type == 'binary') {
    p_snp <- model_summary['snp', 'Pr(>|z|)']
  }

  return(list('beta' = beta_snp,
              'var' = var_snp,
              'p' = p_snp))
}


### acat test for final p-value
acat <- function(p) {
  ### test statistic
  t_acat <- mean(tan((0.5 - p) * pi))
  p_acat <- 0.5 - 1/pi * atan(t_acat)  ### equivalent to pcauchy(t_acat, lower.tail = F)
  return(p_acat)
}

min_p <- function(z, q) {
  ### 2 sided test for P(|Z| >= z)
  t <- max(abs(z))^2
  # p_min <- 1 - pchisq(t, df = 1)^q [equivalent]
  p_min <- 1 - exp(q * pchisq(t, df = 1, log.p = T))
  return(p_min)
}


### Oracle test
### Score test assuming we knew which phenotypes were non-null a priori
oracle_test <- function(score, V_inv, ix_non_null) {
  ### If no non-null pheno-types return p-value of 1
  if(!any(ix_non_null)) {
    return(1)
  }

  ### Filter Score and inverse variance matrix to non-null phenotypes
  score <- score[ix_non_null]
  V_inv <- V_inv[ix_non_null, ix_non_null]

  ### Score test statistic
  q <- sum(ix_non_null)
  t <- as.vector(score %*% V_inv %*% t(t(score)))
  p_score <- pchisq(t, df = q, lower.tail = F)
  return(p_score)
}

oracle_test_uncorrelated <- function(z_score, ix_non_null) {
  ### Score test statistic
  q <- sum(ix_non_null)
  t <- sum(z_score[ix_non_null]^2)
  p_score <- pchisq(t, df = q, lower.tail = F)
  return(p_score)
}

### minp and acat tests for compariong whether one method is better via
### running the test before or after rotation
rotation_tests <- function(score, V, z, q) {
  z_rotated <- z
  z_unrotated <- score/sqrt(diag(V))

  ### min p tests
  minp_rotated <- min_p(z_rotated, q)
  minp_unrotated <- min_p(z_unrotated, q)
  minp_unrotated_bonf <- min(2 * pnorm(-abs(z_unrotated)))

  ### Acat
  acat_rotated <- acat(2 * pnorm(-abs(z_rotated)))
  acat_unrotated <- acat(2 * pnorm(-abs(z_unrotated)))

  df_p <- tibble('minp_rotated' = minp_rotated,
                 'minp_unrotated' = minp_unrotated,
                 'minp_unrotated_bonf' = minp_unrotated_bonf,
                 'acat_rotated' = acat_rotated,
                 'acat_unrotated' = acat_unrotated)

  return(df_p)
}

### ASSET Subset test (based on rotaed Z-scores from mixWAS)
asset_score <- function(component, dataset) {
  ### Get the breakdown of cases + controls
  df_case <-
    map_dfr(dataset$phenotypes, ~{
      as_tibble(.x) %>%
        pivot_longer(cols = everything(),
                     names_to = 'phenotype',
                     values_to = 'case')
    }) %>%
    group_by(phenotype) %>%
    summarise('n_case' = sum(case == 1, na.rm = T),
              'n_control' = sum(case == 0, na.rm = T)) %>%
    arrange(as.numeric(gsub('Y', '', phenotype)))

  results <-
    h.traits(snp.vars = 'SNP_1',
             traits.lab = paste0('Y', 1:component$q),
             beta.hat = component$z,
             sigma.hat = rep(1, component$q),
             side = 2,
             ncase = df_case$n_case,
             ncntl = df_case$n_control,
             meta = F)

  p_value <- h.summary(results)$Subset.2sided$Pvalue


  return(p_value)
}



build_site_cor_mat <- function(phenotypes) {
  ### site specific matrices
  q <- ncol(phenotypes)
  n00 <- matrix(0, nrow = q, ncol = q)
  n10 <- matrix(0, nrow = q, ncol = q)
  n11 <- matrix(0, nrow = q, ncol = q)
  rownames(n00) <- colnames(phenotypes)
  rownames(n10) <- colnames(phenotypes)
  rownames(n11) <- colnames(phenotypes)
  colnames(n00) <- colnames(phenotypes)
  colnames(n10) <- colnames(phenotypes)
  colnames(n11) <- colnames(phenotypes)

  ### Loop over pairs of phenotypes and compute sample overlaps
  for(p1 in colnames(phenotypes)) {
    for(p2 in colnames(phenotypes)) {
      n00[p1, p2] <- sum(phenotypes[,p1] == 0 & phenotypes[,p2] == 0, na.rm = T)
      n10[p1, p2] <- sum(phenotypes[,p1] == 1 & phenotypes[,p2] == 0, na.rm = T)
      n11[p1, p2] <- sum(phenotypes[,p1] == 1 & phenotypes[,p2] == 1, na.rm = T)
    }
  }

  return(list('n00' = n00,
              'n10' = n10,
              'n11' = n11))

}

asset <- function(phenotypes, pheWAS_meta, q) {

  ### Build the correlation matrices needed to run ASSET
  ### N00: Control/Control
  ### N10: Case/Control
  ### N11: Case/Case
  ### Global matrices
  N00 <- matrix(0, nrow = q, ncol = q)
  N10 <- matrix(0, nrow = q, ncol = q)
  N11 <- matrix(0, nrow = q, ncol = q)
  rownames(N00) <- pheWAS_meta$phenotype
  rownames(N10) <- pheWAS_meta$phenotype
  rownames(N11) <- pheWAS_meta$phenotype
  colnames(N00) <- pheWAS_meta$phenotype
  colnames(N10) <- pheWAS_meta$phenotype
  colnames(N11) <- pheWAS_meta$phenotype

  ### Correlation matrices for each individual site
  cor_mats <- map(phenotypes, build_site_cor_mat)

  ### Add components into global correlation matrices
  for(i in 1:length(cor_mats)) {
    ix <- colnames(cor_mats[[i]]$n00)
    N00[ix, ix] <- N00[ix, ix] + cor_mats[[i]]$n00
    N10[ix, ix] <- N10[ix, ix] + cor_mats[[i]]$n10
    N11[ix, ix] <- N11[ix, ix] + cor_mats[[i]]$n11
  }

  ### Run ASSET
  results <-
    h.traits(snp.vars = 'SNP_1',
             traits.lab = pheWAS_meta$phenotype,
             beta.hat = pheWAS_meta$beta_meta,
             sigma.hat = sqrt(pheWAS_meta$var_meta),
             side = 2,
             ncase = diag(N11),
             ncntl = diag(N00),
             cor = list('N00' = N00, 'N10' = N10, 'N11' = N11),
             meta = F)

  p_value <- h.summary(results)$Subset.2sided$Pvalue

  return(p_value)
}


#### Function to modify phenotypes to get healthy controls
get_healthy_controls <- function(phenotypes) {
  ### Get index of patients who are not healthy controls
  ix_not_healthy <- map(phenotypes, ~which(apply(.x, 1, function(vec) {!all(vec == 0)} )))

  ### Make NA to phenotypes
  phenotypes <-
    map2(phenotypes, ix_not_healthy, ~{
      ix_na <- .x[.y,] == 0
      .x[.y,][ix_na] <- NA
      .x
    })

  return(phenotypes)
}

### Wrapper for MultiPhen algorithm
multi_phen <- function(dataset) {
  genotype_matrix <-
    map(dataset$snps, ~{
      G <- matrix(.x) - min(.x)
      rownames(G) <- 1:nrow(G)
      colnames(G) <- 'SNP'
      G
    })

  phenotype_matrix <-
    map2(dataset$phenotypes, dataset$covariates, ~{
      P <- cbind(.x, .y)
      rownames(P) <- 1:nrow(P)
      P
    })

  ### P-Value for Joint Model at Each Site
  pvalues <-
    map2_dbl(genotype_matrix, phenotype_matrix, ~{
      output <-
        MultiPhen::mPhen(genoData = .x[,,drop = F] ,
                         phenoData = .y,
                         covariates = c('PC1', 'PC2', 'PC3', 'PC4', 'age_centered', 'gender'))
      output$Results[1,,,]['JointModel','pvalue']
    })

  ### Deal w/ federation by min P
  return(min(pvalues))
}
