score_test <- function(score, V_inv, q) {
  ### Score test statistic
  t1 <- as.vector(score %*% V_inv %*% t(t(score)))
  p1 <- pchisq(t1, df = q, lower.tail = F)
  return(p1)
}

min_p <- function(score, V_inv, q) {
  V_inv_sqrt <- expm::sqrtm(V_inv)
  z <- V_inv_sqrt %*% score

  ### 2 sided test for P(|Z| >= z)
  p2 <- 2 * pnorm(abs(z), mean = 0, sd = 1, lower.tail = F)
  return(p2)
}

### acat test for final p-value
acat <- function(p) {
  ### test statistic
  t_acat <- mean(tan((0.5 - p) * pi))
  p_acat <- 0.5 - 1/pi * atan(t_acat)  ### equivalent to pcauchy(t_acat, lower.tail = F)
  return(p_acat)
}
