r2 <- read.csv("results/mofacell/r2_total~combined_sn10xLRCluster2.csv")

all_factors_r2 <- r2 %>% dplyr::filter(stringr::str_detect(View, '&', negate = T), Factor == "All")
max(all_factors_r2$R2)
#74.74664
min(all_factors_r2$R2)
#63.30244
r2_batch <- r2 %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2", "Factor7")) %>%
  dplyr::filter(stringr::str_detect(View, "&", negate = TRUE)) %>%
  dplyr::group_by(View) %>%
  dplyr::summarise(R2_Batch = sum(R2))

r2_technical <- r2 %>%
  dplyr::filter(Factor %in% c("Factor6", "Factor4")) %>%
  dplyr::filter(stringr::str_detect(View, "&", negate = TRUE)) %>%
  dplyr::group_by(View) %>%
  dplyr::summarise(R2_Technical = sum(R2))

r2_biological <- r2 %>%
  dplyr::filter(Factor %in% c("Factor3", "Factor5")) %>%
  dplyr::filter(stringr::str_detect(View, "&", negate = TRUE)) %>%
  dplyr::group_by(View) %>%
  dplyr::summarise(R2_Biological = sum(R2))

max(r2_biological$R2_Biological)
# 17.70718
min(r2_biological$R2_Biological)
# 7.784249
r2_outlier <- r2 %>%
  dplyr::filter(Factor %in% stringr::str_c("Factor", 8:20)) %>%
  dplyr::filter(stringr::str_detect(View, "&", negate = TRUE)) %>%
  dplyr::group_by(View) %>%
  dplyr::summarise(R2_Outlier = sum(R2))

df <- r2_batch %>% 
  dplyr::full_join(r2_technical) %>%
  dplyr::full_join(r2_biological) %>%
  dplyr::full_join(r2_outlier) %>%
  tidyr::pivot_longer(cols = starts_with("R2"), names_to = "Type",
                      values_to = 'Sum_R2') %>%
  dplyr::mutate(View = forcats::fct_rev(View)) %>%
  dplyr::mutate(Type = stringr::str_replace(Type, "^R2_", "")) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Biological", "Technical", "Batch", "Outlier")))

pC <- df %>%
  ggplot(aes(x = View, y = Sum_R2)) +
  geom_col(aes(fill = forcats::fct_rev(Type))) + 
  theme_minimal() + 
  coord_flip() + theme_classic() + labs(y = 'R2', fill = "Source of variance") + scale_fill_manual(values = c("Batch" = "#FB8072",
                                                                                                              "Technical" = "#FFFFB3",
                                                                                                              "Outlier" = "#BEBADA",
                                                                                                              "Biological" = "#8DD3C7"
                                                                                                              ))

pC
# 
# design <- "
# 111122"
# (
#   pB + pC + plot_layout(design = design)
# )
