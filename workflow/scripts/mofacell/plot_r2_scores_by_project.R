library(magrittr)
library(ggplot2)

if (exists("snakemake")){
  scores <- read.csv(snakemake@input$scores, row.names = 1)
  r2 <- read.csv(snakemake@input$r2)
  r2_total <- read.csv(snakemake@input$r2_total)
  metadata <- read.csv(snakemake@input$metadata, row.names = 1)
  plt <- snakemake@output$plt
  ggs <- snakemake@output$ggs
} else {
  scores <- read.csv('results/mofacell/scores~Julio_OCEAN_Nereid~sn10x.csv', row.names = 1)
  r2 <- read.csv("results/mofacell/r2_by_sample~Julio_OCEAN_Nereid~sn10x.csv")
  metadata <- read.csv('data/Metadata_For_Charlotte.csv', row.names = 1)
  plt <- 'plots/mofacell/scores_r2_boxplot_by_project~Julio_OCEAN_Nereid~sn10x.pdf'
  ggs <- 'results/mofacell/scores_r2_boxplot_by_project~Julio_OCEAN_Nereid~sn10x.rds'
}

metadata <- metadata %>%
  tibble::rownames_to_column("Sample")

r2 <- r2 %>% dplyr::rename(Sample = Group) %>%
  dplyr::filter(!Factor == 'All')
to_plot <- r2 %>%
  dplyr::left_join(scores, by = c('Sample', 'Factor')) %>%
  dplyr::left_join(metadata, by = c("Sample"))

to_plot_classic <- to_plot %>%
  dplyr::filter(!stringr::str_detect(View, "&"))

views <- unique(to_plot_classic$View)

to_plot_liana <- lapply(views, function(view){
  df <- to_plot %>%
    dplyr::filter(stringr::str_detect(View,
                                      stringr::str_c(view, "&")))
  if (nrow(df) > 0){
    return(df)
  }
  return(NULL)
  }) %>% setNames(views)
# Remove NULL elements
to_plot_liana <- to_plot_liana[lengths(to_plot_liana) != 0]


plot_factor_r2 <- function(f, to_plot){
  # Find best explained view by factor across all samples
  bev <- r2_total %>% dplyr::filter(Factor == f) %>%
    dplyr::slice_max(order_by = R2, n = 1) %>%
    dplyr::select(View) %>% unlist()
  # Drop samples who have is.na(R2) for the best explained view(s)
  drop_samples <- to_plot %>%
    dplyr::filter(Factor == f) %>%
    dplyr::filter((View %in% bev) & (is.na(R2))) %>%
    dplyr::select(Sample) %>% dplyr::distinct() %>% unlist()
  p <- to_plot %>%
    dplyr::filter(!Sample %in% drop_samples) %>%
    dplyr::filter(Factor == f) %>%
    ggplot(aes(x = Project1, y = Score)) + 
    facet_wrap(vars(View), nrow = 3) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = R2), width = 0.3, height = 0, size = 1) +
    scale_color_viridis_c() +
    coord_flip() +
    theme_minimal() +
    labs(y = "Factor score", 
         x = "Project",
         subtitle = f)
  return(p)
}

my_plots_classic <- lapply(unique(to_plot$Factor), plot_factor_r2, to_plot_classic)

if (length(to_plot_liana) > 0){
  # Separate LIANA plots by Factor and by sender cell
  my_plots_liana <- lapply(unique(to_plot$Factor),
                           function(f){
                                        lapply(to_plot_liana,
                                               function(df){
                                                            plot_factor_r2(f, df)})
                                        })
  # Sort all plots by Factor
  my_plots <- lapply(1:length(my_plots_classic), function(i){
    c(my_plots_classic[i], my_plots_liana[[i]])
  }) %>% purrr::flatten()
} else {
  my_plots <- my_plots_classic
}


pdf(plt, width = 12)
my_plots
dev.off()

saveRDS(my_plots, file = ggs)
