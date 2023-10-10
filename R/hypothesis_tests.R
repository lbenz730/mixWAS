score_test <- function(score, V_inv, q) {
  ### Score test statistic
  t <- as.vector(score %*% V_inv %*% t(t(score)))
  p_score <- pchisq(t, df = q, lower.tail = F)
  return(p_score)
}

compute_z_scores <- function(score, V_inv) {
  V_inv_sqrt <- expm::sqrtm(V_inv)
  z <- V_inv_sqrt %*% score
  return(z)
}

min_p <- function(z, q) {
  ### 2 sided test for P(|Z| >= z)
  t <- max(abs(z))^2
  # p_min <- 1 - pchisq(t, df = 1)^q [equivalent]
  p_min <- 1 - exp(q * pchisq(t, df = 1, log.p = T))
  return(p_min)
}

### acat test for final p-value
acat <- function(p) {
  ### test statistic
  t_acat <- mean(tan((0.5 - p) * pi))
  p_acat <- 0.5 - 1/pi * atan(t_acat)  ### equivalent to pcauchy(t_acat, lower.tail = F)
  return(p_acat)
}

#' run hypothesis tests from mixWAS components
#' @param score score vector from mixWAS intermediate output
#' @param V_inv Inverse Varariance Matrix from mixWAS intermediate output
#' @param z Standardized Z-scores from mixWAS intermediate output
#' @param q # of phenotypes
#' @return p-values of aggregate tests
#' @export
run_hypothesis_tests <- function(score, V_inv, z, q) {

  p_score <- score_test(score, V_inv, q) ### Score Test
  p_acat <- acat(2 * pnorm(-abs(z))) ### ACAT of p-values
  p_snp <- acat(c(p_score, p_acat))

  ### Compile List
  p_values <-
    list('p_score' = p_score,
         'p_acat' = p_acat,
         'p_snp' = p_snp)

  return(p_values)
}

