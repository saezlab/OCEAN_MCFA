library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(Hmisc)

scores <- read.csv('results/robustness/scores~Julio_OCEAN_Nereid~sn10xcanon2.csv')[,-c(1)]
loadings <- read.csv('results/robustness/loadings~Julio_OCEAN_Nereid~sn10xcanon2.csv')[,-c(1)]
r2_sample <- read.csv("results/robustness/r2_by_sample~Julio_OCEAN_Nereid~sn10xcanon2.csv") %>%
     dplyr::rename('Sample' = 'Group')
r2_total <- read.csv("results/robustness/r2_total~Julio_OCEAN_Nereid~sn10xcanon2.csv")

# Find best explained view by model and factor across all samples
bev <- r2_total %>% dplyr::group_by(n_factors, Factor) %>%
  dplyr::slice_max(order_by = R2, n = 1) %>%
  dplyr::select(n_factors, Factor, View)
# Drop samples who have is.na(R2) for the best explained view(s)
drop_samples <- r2_sample %>%
  dplyr::semi_join(bev, by = c("n_factors", "Factor", "View")) %>%
  dplyr::filter(is.na(R2)) %>%
  dplyr::select(n_factors, Factor, Sample, View) %>%
  dplyr::distinct()

scores_filtered <- scores %>%
  dplyr::anti_join(drop_samples, by = c("n_factors", "Factor", "Sample"))

#####
# Factor robustness
#####
pivoted_scores_all <- scores_filtered %>%
  tidyr::pivot_wider(names_from = c(n_factors,
                                    Factor),
                     values_from = Score,
                     values_fill = NULL)

pivoted_loadings_all <- loadings %>%
  tidyr::pivot_wider(names_from = c(n_factors, 
                                    Factor),
                     values_from = Loading,
                     values_fill = NULL)

mat_scores_all <- pivoted_scores_all %>%
  dplyr::select(-c(Sample)) %>%
  as.matrix()

corr_scores_all <- Hmisc::rcorr(mat_scores_all,
                                type="pearson")$r

mat_loadings_all <- pivoted_loadings_all %>%
  dplyr::select(-c(view, variable)) %>%
  as.matrix()

corr_loadings_all <- Hmisc::rcorr(mat_loadings_all,
                                  type="pearson")$r


all_correlations <- corr_scores_all*corr_loadings_all

df <- all_correlations %>%
  data.frame %>%
  rownames_to_column("factor1") %>%
  dplyr::mutate(factor1 = stringr::str_c("X", factor1)) %>%
  tidyr::separate(col = factor1, into = c("model1", "factor1")) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = c("model2", "factor2"),
               names_sep = "_",
               values_to = "correlation")

df2 <- df %>%
  dplyr::distinct() %>%
  dplyr::filter(!(model1 == model2 & factor1 == factor2)) %>%
  dplyr::group_by(factor1, model1, model2) %>%
  dplyr::summarise(max_corr = max(correlation)) %>%
  ungroup() %>%
  dplyr::mutate(model1 = as.integer(str_extract(model1, "[:digit:]+")),
                model2 = as.integer(str_extract(model2, "[:digit:]+")))

# Take the mean of the max correlation amongst all models which are bigger:
df3 <- df2 %>%
  dplyr::group_by(factor1, model1) %>%
  dplyr::filter(model2 > model1) %>%
  dplyr::summarise(robustness = mean(max_corr))

# Take the mean across the max correlation of all other models:  
factor_robustness <- df2 %>%
  dplyr::group_by(factor1, model1) %>%
  dplyr::summarise(robustness = mean(max_corr)) %>%
  dplyr::rename(factor = factor1, model = model1)

factor_mean_r2 <- r2_sample %>%
  dplyr::rename(factor = Factor, model = n_factors) %>%
  dplyr::filter(factor != "All") %>%
  dplyr::group_by(factor, model) %>%
  dplyr::summarise(mean_R2 = mean(R2, na.rm = T))

mat <- all_correlations
HM <- Heatmap(mat,
              row_km = 15,
              row_km_repeats = 500)#, column_km = 20, column_km_repeats = 50)
HM <- draw(HM)
clustering <- row_order(HM)
factor_clusters <- lapply(1:length(clustering),
                          function(i) row.names(mat[clustering[[i]],])) %>%
  setNames(1:length(clustering)) %>%
  map2_df(seq_along(.), ~ tibble(factor = .x, cluster = .y)) %>%
  tidyr::separate(col = factor, into = c("model", "factor")) %>%
  dplyr::mutate(model = as.integer(model))

df <- purrr::reduce(list(factor_clusters, factor_mean_r2, factor_robustness),
                    dplyr::inner_join,
                    by = c("model", "factor")) 

p1 <- df %>%
  ggplot(aes(x = mean_R2, y = robustness, colour = factor(model))) +
  geom_point() + theme_minimal() + scale_color_discrete("Model")

p2 <- df %>% ggplot(aes(x = mean_R2, y = robustness, colour = factor(cluster))) +
  geom_point() + theme_minimal() + scale_color_discrete("Cluster")

p1 + p2
#####
# Correlation between factor scores
#####


pivoted_scores <- scores %>%
  dplyr::filter(n_factors >=10 & n_factors <=15) %>%
  tidyr::pivot_wider(names_from = c(n_factors, Factor),
                     values_from = Score,
                     values_fill = NULL)

mat_scores <- pivoted_scores %>%
  dplyr::select(-c(Sample)) %>%
  as.matrix()

corr_scores <- Hmisc::rcorr(mat_scores, type="pearson")$r

#####
# Correlation between factor loadings
#####


pivoted_loadings <- loadings %>%
  dplyr::filter(n_factors >=10 & n_factors <=15) %>%
  tidyr::pivot_wider(names_from = c(n_factors, Factor),
                     values_from = Loading,
                     values_fill = NULL)

mat_loadings <- pivoted_loadings %>%
  dplyr::select(-c(view, variable)) %>%
  as.matrix()

corr_loadings <- Hmisc::rcorr(mat_loadings, type="pearson")$r


####
# Combination of scores and loadings
####
mat <- corr_scores*corr_loadings

add_stars_make_diag <- function(j, i, x, y, w, h, fill){
  if(i <= j) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    if (mat[i, j] >= 0.95^2) {
      grid.text("*", x, y)
    }
  }
}

add_stars <- function(j, i, x, y, w, h, fill){
  if (mat[i, j] >= 0.95^2) {
    grid.text("*", x, y)
  }
}

f1 <- circlize::colorRamp2(seq(-1, 1, length = 3),
                           c("blue", "#EEEEEE", "red"))

colnames(mat) <- NULL

Heatmap(mat,
        rect_gp = gpar(type = "none"), 
        cell_fun = add_stars_make_diag,
        row_names_gp = gpar(fontsize = 8),
        #column_names_gp = gpar(fontsize = 8),
        cluster_rows = F,
        cluster_row_slices = F,
        #cluster_columns = F,
        cluster_column_slices = F,
        row_split = rep(10:15, times = 10:15),
        row_title = NULL,
        column_split = rep(10:15, times = 10:15),
        row_names_side = "right",
        cluster_columns = F,
        col = f1,
        heatmap_legend_param = list(title = "r"))


####
# R2
####


r2_all <- r2_sample %>% dplyr::filter(Factor == "All")

n_sample_views <- r2_all %>%
  tidyr::drop_na() %>%
  dplyr::select(Sample, View) %>%
  dplyr::distinct() %>% nrow()

thresholds <- c(0, 25, 50)

to_plot <- lapply(thresholds, function(threshold){
  df <- r2_all %>% 
    dplyr::filter(R2 >= threshold) %>%
    dplyr::group_by(n_factors) %>%
    dplyr::summarize(cumul_prop = n_distinct(paste(View,
                                                   Sample,
                                                   sep = "_")) / n_sample_views)
  df$threshold <- threshold
  return(df)
})

to_plot <- dplyr::bind_rows(to_plot)

ggplot(to_plot, aes(x = n_factors, y = cumul_prop, color = factor(threshold))) +
  geom_line() +
  labs(x = "Number of factors in model",
       y = "% views reaching R2 threshold",
       colour = "R2 threshold") +
  ylim(0,1) +  theme_minimal()

df <- r2_all %>%
  dplyr::group_by(View, Factor, n_factors) %>%
  dplyr::summarise(median_R2 = median(R2, na.rm = T),
                   mean_R2 = mean(R2, na.rm = T))

library(ggrepel)
df %>% dplyr::filter(Factor == 'All') %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(View),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = median_R2,
             color = View)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 3.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,22.5) +
  labs(x = "Number of factors in model",
       y = expression("Median"~ R^2 ~ "across samples"))

df %>% dplyr::filter(Factor == 'All') %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(View),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = mean_R2,
             color = View)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 3.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,22.5) +
  labs(x = "Number of factors in model",
       y = expression("Mean"~ R^2 ~ "across samples"))
