library(ggplot2)
library(umap)
library(magrittr)
library(dplyr)
library(patchwork)

if (exists("snakemake")){
  scores <- read.csv(snakemake@input$scores, row.names = 1)
  metadata <- read.csv(snakemake@input$metadata, row.names = 1)
  data <- snakemake@output$data
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else {
  scores <- read.csv('results/mofacell/scores~Julio_OCEAN_Nereid~all.csv', row.names = 1)
  metadata <- read.csv('data/Metadata_For_Charlotte.csv', row.names = 1)
  data <- "results/mofacell/umap~Julio_OCEAN_Nereid~all.csv"
  plt <- "plots/mofacell/umap~Julio_OCEAN_Nereid~all.pdf"
  gg <- "results/mofacell/umap~Julio_OCEAN_Nereid~all.Rds"
}

metadata <- metadata %>%
  tibble::rownames_to_column("Sample")

to_plot <- scores %>%
  #dplyr::filter(n_factors == 4) %>%
  #dplyr::select(-c(n_factors)) %>%
  tidyr::pivot_wider(names_from = Factor, values_from = Score) %>%
  tibble::column_to_rownames("Sample") %>% as.matrix

# Run UMAP
umap_result <- umap(to_plot, n_neighbors = 25, min_dist = 0.5, n_threads = 1)

# Access the UMAP coordinates
umap_coordinates <- umap_result$layout

simplify_cohort <- function(cohort){
  ifelse(stringr::str_detect(cohort, "IgAN"), "IgAN",
         ifelse(cohort %in% c("AIN/T2D", "C1QN", "C3_Glomerulopathy",
                              "COVAN", "Immune_Cmpx_GN",
                              "MN_Vax", "MPGN"), "Other",
                ifelse(stringr::str_detect(cohort, "FSGS"), "FSGS", cohort)))
}
# Create a data frame with UMAP coordinates
umap_df <- as.data.frame(umap_coordinates) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  tibble::rownames_to_column("Sample") %>% dplyr::left_join(metadata) %>%
  dplyr::mutate(cohort_simplified = simplify_cohort(Cohort))

# # Plot the UMAP result
# (
#   p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = Cohort)) +
#     geom_point(size = 2) +
#     labs(colour = "Disease group") +
#     #scale_colour_brewer(palette = "Set2") + 
#     theme_minimal()
# )

(
  p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cohort_simplified)) +
    geom_point(size = 2) +
    labs(subtitle = "Disease", colour = "Disease group") +
    scale_colour_brewer(palette = "Set2") + 
    theme_minimal()
)

(
  p1a <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = Prep)) +
    geom_point(size = 2) +
    labs(subtitle = "Omic", colour = "Omic") +
    scale_colour_brewer(palette = "Set2") + 
    theme_minimal()
)

(
  p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = as.numeric(eGFRatBx_NEPTUNE))) +
    geom_point(size = 2) +
    labs(subtitle = "eGFR", colour = expression(paste("eGFR at biopsy\n(ml/min/1.73", m^2, ")"))) +
    scale_colour_viridis_c() + theme_minimal()
)

(
  p3 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = as.numeric(Age))) +
    geom_point(size = 2) +
    labs(colour = "Age", subtitle = "Age") +
    scale_colour_viridis_c() + theme_minimal()
)

(
  p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = log10(as.numeric(UPCRatBx_NEPTUNE)))) +
    geom_point(size = 2) +
    labs(colour = "log10(uPCR)", subtitle = "Proteinuria") +
    scale_colour_distiller(palette = "Spectral", limits = c(-2,2)) + theme_minimal()
)

(
  p <- p1 + p1a + p2 + p4 + 
  plot_layout(guides = 'collect') #&
)
  # theme(legend.position = 'bottom',
  #       legend.direction = 'vertical')

ggsave(plot = p, filename = plt, device = "pdf", width = 30, height = 22, units = "cm")
saveRDS(p, file = gg)
write.csv(umap_df, file = data)
