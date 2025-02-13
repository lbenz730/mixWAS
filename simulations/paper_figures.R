library(tidyverse)
library(patchwork)

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 12),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


### Updated Main Figures
### Final 1: Same Signs
df_sims <-
  bind_rows(
    read_csv('outputs/v3/simulation_run_v3_same_sign_1/combined.csv') %>% mutate('sigma' = 'Positive Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_same_sign_2/combined.csv') %>% mutate('sigma' = 'No Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_same_sign_3/combined.csv') %>% mutate('sigma' = 'Negative Phenotype Correlation')
  )

df_final_1 <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin + n_true_con, 'Non-Null Phenotypes'), n_true_bin + n_true_con)) %>%
  mutate('direction' = gsub('\\s+--','', paste(direction_bin, direction_con))) %>%
  mutate('direction' = case_when(direction == 'Same Positive Same Positive' ~ 'Same Direction',
                                 direction == 'Same Positive Opposite' ~ 'Opposite Direction',
                                 direction == 'Same Positive Same Negative' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction,
                                   'Same Direction',
                                   'Opposite Direction'),
         'sigma' = fct_relevel(sigma,
                               'Positive Phenotype Correlation',
                               'No Phenotype Correlation',
                               'Negative Phenotype Correlation')) %>%
  group_by(method, power, sigma, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen', 'Oracle'))

ggplot(df_final_1, aes(x = max_beta_bin, y = power)) +
  facet_grid(sigma~sparsity, labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(0, 0.2)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Mixed Datatype Phenotypes | Same Direction Effects (Positive) | MAF: 20% | Prevelance: 30%',
       color = '')

ggsave('paper_figures/final_simulation_same.png', width = 16/1.2, height = 9/1.2)


### Final 2: Opp Signs
df_sims <-
  bind_rows(
    read_csv('outputs/v3/simulation_run_v3_opp_sign_1/combined.csv') %>% mutate('sigma' = 'Positive Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_opp_sign_2/combined.csv') %>% mutate('sigma' = 'No Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_opp_sign_3/combined.csv') %>% mutate('sigma' = 'Negative Phenotype Correlation')
  )

df_final_2 <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin + n_true_con, 'Non-Null Phenotypes'), n_true_bin + n_true_con)) %>%
  mutate('direction' = gsub('\\s+--','', paste(direction_bin, direction_con))) %>%
  mutate('direction' = case_when(direction == 'Same Positive Same Positive' ~ 'Same Direction',
                                 direction == 'Same Positive Opposite' ~ 'Opposite Direction',
                                 direction == 'Same Positive Same Negative' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction,
                                   'Same Direction',
                                   'Opposite Direction'),
         'sigma' = fct_relevel(sigma,
                               'Positive Phenotype Correlation',
                               'No Phenotype Correlation',
                               'Negative Phenotype Correlation')) %>%
  group_by(method, power, sigma, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen', 'Oracle'))

ggplot(df_final_2, aes(x = max_beta_bin, y = power)) +
  facet_grid(sigma~sparsity, labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  # scale_color_manual(values = gg_color_hue(4)[c(1:2, 4)]) +
  scale_x_continuous(limits = c(0, 0.2)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Mixed Datatype Phenotypes | Opposite Direction Effects (Continuous Phenotypes) | MAF: 20% | Prevelance: 30%',
       color = '')

ggsave('paper_figures/final_simulation_opp.png', width = 16/1.2, height = 9/1.2)


### Correlation Matrix
same_sim_inputs1 <- read_rds('inputs/v3/simulation_run_v3_same_sign_1.rds')
same_sim_inputs2 <- read_rds('inputs/v3/simulation_run_v3_same_sign_2.rds')
same_sim_inputs3 <- read_rds('inputs/v3/simulation_run_v3_same_sign_3.rds')
same_params1 <- same_sim_inputs1$params
same_params2 <- same_sim_inputs2$params
same_params3 <- same_sim_inputs3$params

opp_sim_inputs1 <- read_rds('inputs/v3/simulation_run_v3_opp_sign_1.rds')
opp_sim_inputs2 <- read_rds('inputs/v3/simulation_run_v3_opp_sign_2.rds')
opp_sim_inputs3 <- read_rds('inputs/v3/simulation_run_v3_opp_sign_3.rds')
opp_params1 <- opp_sim_inputs1$params
opp_params2 <- opp_sim_inputs2$params
opp_params3 <- opp_sim_inputs3$params

bind_rows(
  as_tibble(same_params1$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 1: Same Direction Effects',
           'sigma' = 'Positive Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation'),
  as_tibble(same_params2$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 1: Same Direction Effects',
           'sigma' = 'No Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation'),
  as_tibble(same_params3$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 1: Same Direction Effects',
           'sigma' = 'Negative Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation'),
  as_tibble(opp_params1$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 2: Opposite Direction Effects',
           'sigma' = 'Positive Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation'),
  as_tibble(opp_params2$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 2: Opposite Direction Effects',
           'sigma' = 'No Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation'),
  as_tibble(opp_params3$Sigma) %>%
    set_names(paste0('Y', 1:8)) %>%
    mutate('phenotype' = paste0('Y', 1:8),
           'simulation' = 'Simulation 2: Opposite Direction Effects',
           'sigma' = 'Negative Phenotype Correlation') %>%
    pivot_longer(cols = starts_with('Y'),
                 names_to = 'phenotype2',
                 values_to = 'correlation')
) %>%
  mutate('sigma' = fct_relevel(sigma,
                               'Positive Phenotype Correlation',
                               'No Phenotype Correlation',
                               'Negative Phenotype Correlation')) %>%
  ggplot(aes(x = phenotype, y = fct_reorder(phenotype2, rep(8:1, 8 * 6)))) +
  facet_grid(simulation~sigma, labeller = label_wrap_gen(width = 25)) +
  geom_tile(aes(fill = as.factor(correlation)), col = 'black',  alpha = 0.8) +
  scale_fill_manual(values = c('red', 'white', 'dodgerblue', 'black')) +
  theme(panel.grid.major = element_blank()) +
  labs(x = '',
       y = '',
       fill = 'Correlation')

ggsave('paper_figures/correlation_simulation_main.png', height = 6, width = 9)

### Appendix figures
### Simulation 2
sim_inputs <- read_rds('inputs/v3/simulation_run_v3_2.rds')
params <- sim_inputs$params
df_sims <- read_csv('outputs/v3/simulation_run_v3_2/combined.csv')


df2 <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin, 'Non-Null Phenotypes'), n_true_bin)) %>%
  mutate('direction' = gsub('\\s+--','', direction_bin)) %>%
  mutate('direction' = case_when(direction == 'Same Positive' ~ 'Same Direction (Positive)',
                                 direction == 'Same Negative' ~ 'Same Direction (Negative)',
                                 direction == 'Opposite' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction, 'Same Direction (Positive)', 'Same Direction (Negative)', 'Opposite Direction')) %>%
  group_by(method, power, direction, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen', 'ASSET (Subset Search)', 'Oracle'))

ggplot(df2, aes(x = max_beta_bin, y = power)) +
  facet_grid(direction~sparsity, labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(0, 0.2)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Binary Phenotypes | Varying Sparsity and Direction | MAF: 20% | Prevalence: 30%',
       color = '')

ggsave('paper_figures/simulation_2.png', width = 16/1.25, height = 9/1.25)


### Simulation 3
sim_inputs <- read_rds('inputs/v3/simulation_run_v3_3.rds')
params <- sim_inputs$params
df_sims <- read_csv('outputs/v3/simulation_run_v3_3/combined.csv')


df3 <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin, 'Non-Null Phenotypes'), n_true_bin)) %>%
  mutate('direction' = gsub('\\s+--','', direction_bin)) %>%
  mutate('direction' = case_when(direction == 'Same Positive' ~ 'Same Direction (Positive)',
                                 direction == 'Same Negative' ~ 'Same Direction (Negative)',
                                 direction == 'Opposite' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction, 'Same Direction (Positive)', 'Same Direction (Negative)', 'Opposite Direction')) %>%
  group_by(method, power, direction, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen', 'ASSET (Subset Search)', 'Oracle'))

ggplot(df3, aes(x = max_beta_bin, y = power)) +
  facet_grid(direction~sparsity, scales = 'free_x', labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Binary Phenotypes | Varying Sparsity and Direction | MAF: 5% | Prevalence: 10%',
       color = '')

ggsave('paper_figures/simulation_3.png', width = 16/1.25, height = 9/1.25)

as_tibble(params$Sigma) %>%
  set_names(paste0('Y', 1:8)) %>%
  mutate('phenotype' = paste0('Y', 1:8)) %>%
  pivot_longer(cols = starts_with('Y'),
               names_to = 'phenotype2',
               values_to = 'correlation') %>%
  ggplot(aes(x = phenotype, y = fct_reorder(phenotype2, rep(8:1, 8)))) +
  geom_tile(aes(fill = as.factor(correlation)), col = 'black',  alpha = 0.8) +
  scale_fill_manual(values = c('white', 'orange', 'pink', 'black')) +
  theme(panel.grid.major = element_blank()) +
  labs(x = '',
       y = '',
       fill = 'Correlation')
ggsave('paper_figures/correlation_simulation_3.png', width = 16/2, height = 9/2)


### Simulation 4
sim_inputs <- read_rds('inputs/v3/simulation_run_v3_4.rds')
params <- sim_inputs$params
df_sims <- read_csv('outputs/v3/simulation_run_v3_4/combined.csv')


df4 <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin, 'Non-Null Phenotypes'), n_true_bin)) %>%
  mutate('direction' = gsub('\\s+--','', direction_bin)) %>%
  mutate('direction' = case_when(direction == 'Same Positive' ~ 'Same Direction (Positive)',
                                 direction == 'Same Negative' ~ 'Same Direction (Negative)',
                                 direction == 'Opposite' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction, 'Same Direction (Positive)', 'Same Direction (Negative)', 'Opposite Direction')) %>%
  group_by(method, power, direction, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'ASSET (Subset Search)', 'Oracle'))

ggplot(df4, aes(x = max_beta_bin, y = power)) +
  facet_grid(direction~sparsity, scales = 'free_x', labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(0, 0.3)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Binary Phenotypes | Varying Sparsity and Direction | MAF: 20% | Prevalence: 30%\n Healthy Controls',
       color = '')

ggsave('paper_figures/simulation_4.png', width = 16/1.25, height = 9/1.25)

as_tibble(params$Sigma) %>%
  set_names(paste0('Y', 1:8)) %>%
  mutate('phenotype' = paste0('Y', 1:8)) %>%
  pivot_longer(cols = starts_with('Y'),
               names_to = 'phenotype2',
               values_to = 'correlation') %>%
  ggplot(aes(x = phenotype, y = fct_reorder(phenotype2, rep(8:1, 8)))) +
  geom_tile(aes(fill = factor(correlation, levels = c(0, 0.4, 0.7, 1))), col = 'black',  alpha = 0.8) +
  scale_fill_manual(values = c('white', 'orange', 'pink', 'black'), drop = F) +
  theme(panel.grid.major = element_blank()) +
  labs(x = '',
       y = '',
       fill = 'Correlation')
ggsave('paper_figures/correlation_simulation_4.png', width = 16/2, height = 9/2)


### Simulation 5 (MAR)
df_sims <-
  bind_rows(
    read_csv('outputs/v3/simulation_run_v3_mar_1/combined.csv') %>% mutate('sigma' = 'Positive Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_mar_2/combined.csv') %>% mutate('sigma' = 'No Phenotype Correlation'),
    read_csv('outputs/v3/simulation_run_v3_mar_3/combined.csv') %>% mutate('sigma' = 'Negative Phenotype Correlation')
  )

df_final_mar <-
  df_sims %>%
  pivot_longer(cols = starts_with('p_'),
               names_to = 'method',
               values_to = 'power') %>%
  mutate('method' = case_when(method == 'p_snp' ~ 'mixWAS',
                              method == 'p_pheWAS_mega' ~ 'PheWAS Mega',
                              method == 'p_pheWAS_meta' ~ 'PheWAS (min P)',
                              method == 'p_score' ~ 'mixWAS (Score Only)',
                              method == 'p_acat' ~ 'ACAT of Score P-Values',
                              method == 'p_pheWAS_mega_acat' ~ 'PheWAS Mega ACAT',
                              method == 'p_pheWAS_meta_acat' ~ 'PheWAS Meta ACAT',
                              method == 'p_pheWAS_mega_acat_p' ~ 'ACAT of PheWAS Mega P-Values',
                              method == 'p_pheWAS_meta_acat_p' ~ 'PheWAS (ACAT)',
                              method == 'p_oracle_uncorrelated' ~ 'Oracle',
                              method == 'p_asset_meta' ~ 'ASSET (Subset Search)',
                              method == 'p_multiphen' ~ 'MultiPhen',
                              method == 'p_pheWAS_hmp' ~ 'PheWAS (HMP)'
  )) %>%
  mutate('method' = fct_relevel(method, 'mixWAS', 'mixWAS (Score Only)', 'ACAT of Score P-Values', 'PheWAS (min P)', 'PheWAS Mega',
                                'PheWAS Mega ACAT', 'PheWAS Meta ACAT',
                                'ACAT of PheWAS Mega P-Values',  'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen',
                                'ASSET (Subset Search)', 'Oracle')) %>%
  mutate('sparsity' = fct_reorder(paste(n_true_bin + n_true_con, 'Non-Null Phenotypes'), n_true_bin + n_true_con)) %>%
  mutate('direction' = gsub('\\s+--','', paste(direction_bin, direction_con))) %>%
  mutate('direction' = case_when(direction == 'Same Positive Same Positive' ~ 'Same Direction',
                                 direction == 'Same Positive Opposite' ~ 'Opposite Direction',
                                 direction == 'Same Positive Same Negative' ~ 'Opposite Direction')) %>%
  mutate('direction' = fct_relevel(direction,
                                   'Same Direction',
                                   'Opposite Direction'),
         'sigma' = fct_relevel(sigma,
                               'Positive Phenotype Correlation',
                               'No Phenotype Correlation',
                               'Negative Phenotype Correlation')) %>%
  group_by(method, power, sigma, sparsity) %>%
  mutate('min_beta' = min(max_beta_bin)) %>%
  filter(max_beta_bin == min_beta) %>%
  ungroup() %>%
  filter(method %in% c('mixWAS', 'PheWAS (min P)', 'PheWAS (ACAT)', 'PheWAS (HMP)', 'MultiPhen', 'Oracle'))

ggplot(df_final_mar, aes(x = max_beta_bin, y = power)) +
  facet_grid(sigma~sparsity, labeller = label_wrap_gen(width = 20)) +
  geom_line(aes(col = method)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(0, 0.2)) +
  labs(x = expression(paste('Effect Size (', beta, ')')),
       y = 'Power',
       title = 'Power for Cross-Phenotype Association Test',
       subtitle = 'Mixed Datatype Phenotypes | Same Direction Effects (Positive) | MAF: 20% | Prevelance: 30% | MAR Outcomes',
       color = '')

ggsave('paper_figures/simulation_5.png', width = 16/1.2, height = 9/1.2)
