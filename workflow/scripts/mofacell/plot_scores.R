library(magrittr)
library(ggplot2)

if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else {
  data <- read.csv('results/mofacell/scores~Julio_OCEAN_Nereid~all.csv', row.names = 1)
  metadata <- read.csv('data/Metadata_For_Charlotte.csv', row.names = 1)
  plt <- 'plots/mofacell/scores_r2_boxplot~Julio_OCEAN_Nereid~all.pdf'
}

simplify_cohort <- function(cohort){
  ifelse(stringr::str_detect(cohort, "IgAN"), "IgAN",
         ifelse(cohort %in% c("AIN/T2D", "C1QN", "C3_Glomerulopathy",
                              "COVAN", "Immune_Cmpx_GN",
                              "MN_Vax", "MPGN"), "Other",
                ifelse(stringr::str_detect(cohort, "FSGS"), "FSGS", cohort)))
}
metadata <- metadata %>%
  tibble::rownames_to_column("Sample")

to_plot <- data %>%
  dplyr::rename(ID = Sample) %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(cohort_simplified = simplify_cohort(Cohort))

# Factor 1
(
  p <- to_plot %>%
    dplyr::filter(Factor == 'Factor5') %>%
    dplyr::select(Cohort, Score) %>%
    dplyr::distinct() %>%
    ggplot(aes(x = Cohort, y = Score, colour = Cohort)) + 
    #scale_colour_brewer(palette = "Set2") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = Cohort), width = 0.3, height = 0) +
    theme_minimal() +
    labs(y = "Factor score", 
         x = "Disease group") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)

(
  p1 <- to_plot %>%
    dplyr::select(ID, Cohort, Score, Factor) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(id_cols = c(ID, Cohort),
                       names_from = Factor,
                       values_from = Score) %>%
    ggplot(aes(x = Factor1, y = Factor3, colour = cohort_simplified)) +
    geom_point() +
    #scale_colour_brewer(palette = "Set2") 
    theme_minimal() +
    labs(colour = "Disease group") +
    theme(legend.position = 'top', legend.direction = "horizontal")
)



# ggsave(plot = p1, filename = plt1, device = "pdf", width = 15, height = 12, units = "cm")
# ggsave(plot = p2, filename = plt2, device = "pdf", width = 15, height = 12, units = "cm")
# ggsave(plot = p3, filename = plt3, device = "pdf", width = 15, height = 12, units = "cm")
# ggsave(plot = p4, filename = plt4, device = "pdf", width = 15, height = 12, units = "cm")

