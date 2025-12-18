
hallmark_res <- read.csv('results/downstream/hallmark_factor5.csv')[,-1]
# Flip Factor 5 for consistency with preprint (sign invariant)
df <- hallmark_res %>%
  dplyr::select(Term, FDR.p.value, view, positive) %>%
  dplyr::mutate(score = ifelse(positive == 'False', # Normally 'True'
                               -log10(FDR.p.value),
                               log10(FDR.p.value)))
terms_to_plot <- df %>%
  dplyr::filter(FDR.p.value < 1e-5) %>%
  dplyr::select(Term) %>% unlist() %>% unique()

mat <- df %>%
  dplyr::filter(Term %in% terms_to_plot) %>%
  dplyr::transmute(Term = stringr::str_replace(Term, "HALLMARK_", ""), view, score, positive) %>%
  dplyr::group_by(Term, view) %>%                     # Two directions: take the one with 
  dplyr::slice_max(order_by = abs(score), n = 1) %>%  # max absolute score (most relevant)
  dplyr::ungroup() %>%
  dplyr::select(-positive) %>%
  tidyr::pivot_wider(names_from = view, values_from = score, values_fill = 0) %>%
  tibble::column_to_rownames('Term') %>% as.matrix()

kclus <- kmeans(mat, 2)

cn = colnames(mat)
split = factor(paste0("Cluster", kclus$cluster), levels = c("Cluster1", "Cluster2"))
col_fun = circlize::colorRamp2(c(-max(abs(mat)), 0, max(abs(mat))), c("blue", "white", "red"))
orderedHM <- ComplexHeatmap::Heatmap(mat,
                                     #rect_gp = gpar(type = "none"),
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     show_row_dend = FALSE,
                                     split = split,
                                     cluster_row_slices = FALSE,
                                     show_row_names = TRUE,
                                     col = col_fun,#viridis::viridis(100),
                                     name = "-log10(padj)",
                                     show_column_names = FALSE, 
                                     bottom_annotation = HeatmapAnnotation(
                                       text = anno_text(cn, rot = 45, location = unit(0.9, "npc"), just = "right"),
                                       annotation_height = max_text_width(cn)
                                     ))
draw(orderedHM, padding = unit(c(2, 10, 2, 60), "mm"), heatmap_legend_side = "bottom")



hallmark_res3 <- read.csv('results/downstream/hallmark_factor3.csv')[,-1]

df <- hallmark_res3 %>%
  dplyr::select(Term, FDR.p.value, view, positive) %>%
  dplyr::mutate(score = ifelse(positive == 'True', 
                               -log10(FDR.p.value),
                               log10(FDR.p.value)))
terms_to_plot <- df %>%
  dplyr::filter(FDR.p.value < 1e-5) %>%
  dplyr::select(Term) %>% unlist() %>% unique()

mat <- df %>%
  dplyr::filter(Term %in% terms_to_plot) %>%
  dplyr::transmute(Term = stringr::str_replace(Term, "HALLMARK_", ""), view, score, positive) %>%
  dplyr::group_by(Term, view) %>%                     # Two directions: take the one with 
  dplyr::slice_max(order_by = abs(score), n = 1) %>%  # max absolute score (most relevant)
  dplyr::ungroup() %>%
  dplyr::select(-positive) %>%
  tidyr::pivot_wider(names_from = view, values_from = score, values_fill = 0) %>%
  tibble::column_to_rownames('Term') %>% as.matrix()

kclus <- kmeans(mat, 3)

cn = colnames(mat)
split = factor(paste0("Cluster", kclus$cluster), levels = c("Cluster1", "Cluster2", "Cluster3"))
col_fun = circlize::colorRamp2(c(-max(abs(mat)), 0, max(abs(mat))), c("blue", "white", "red"))
orderedHM <- ComplexHeatmap::Heatmap(mat,
                                     #rect_gp = gpar(type = "none"),
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     show_row_dend = FALSE,
                                     #split = split,
                                     cluster_row_slices = FALSE,
                                     show_row_names = TRUE,
                                     col = col_fun,#viridis::viridis(100),
                                     name = "Score",
                                     show_column_names = FALSE, 
                                     bottom_annotation = HeatmapAnnotation(
                                       text = anno_text(cn, rot = 45, location = unit(0.9, "npc"), just = "right"),
                                       annotation_height = max_text_width(cn)
                                     ))
draw(orderedHM, padding = unit(c(2, 10, 2, 60), "mm"), heatmap_legend_side = "bottom")

