library(ggplot2)
library(patchwork)
library(magrittr)
library(ComplexHeatmap)
r2 <- read.csv('results/mofacell/r2_total~combined_sn10xLRCluster2.csv')

to_plot <- r2 %>%
  dplyr::filter(Factor %in% c("Factor3", "Factor5")) %>%
  dplyr::filter(!stringr::str_detect(View, "&"))
(
  pH <- ggplot(to_plot, aes(y = forcats::fct_rev(View), x = Factor, fill = R2))
  + geom_tile()
  + scale_fill_viridis_c(option = "plasma", limits = c(0, 15))  + theme_classic()
  + labs(x = element_blank(), y = 'View') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)))
)


acts <- read.csv("results/downstream/progeny_acts~combined_sn10xLRCluster2.csv") %>% setNames(c("Factor", "Pathway", "Score", "View"))
pvals <- read.csv("results/downstream/progeny_pvals~combined_sn10xLRCluster2.csv")  %>% setNames(c("Factor", "Pathway", "pval", "View"))
progeny <- dplyr::full_join(acts,pvals)

progeny_f5 <- progeny %>% dplyr::filter(Factor == "Factor5") %>%
  dplyr::filter(Pathway %in% c("EGFR", "JAK-STAT", "MAPK", "NFkB", "TGFb", "TNFa")) %>%
  dplyr::mutate(Score = -Score) # Invert Factor5 for consistency with preprint (sign invariant)

(
  pI1 <- progeny_f5 %>%
    dplyr::mutate(View = forcats::fct_rev(View),
                  sig = ifelse(pval < 0.01, "*", "")) %>%
    ggplot(aes(x = Pathway, y = View, fill = Score)) +
    geom_tile() + 
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-8, 8)
    )  + theme_classic() +
    geom_text(aes(label = sig), color = "white") +
    labs(y = element_blank(), x = element_blank(), subtitle = 'Factor 5') +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)))
)
# 3.78x3.16

progeny_f3 <- progeny %>% dplyr::filter(Factor == "Factor3") %>%
  dplyr::filter(Pathway %in% c("EGFR", "JAK-STAT", "MAPK", "NFkB", "TGFb", "TNFa"))

(
  pI2 <- progeny_f3 %>%
    dplyr::mutate(View = forcats::fct_rev(View),
                  sig = ifelse(pval < 0.01, "*", "")) %>%
    ggplot(aes(x = Pathway, y = View, fill = Score)) +
    geom_tile() + 
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-8, 8)
    )  + theme_classic() +
    geom_text(aes(label = sig), color = "white") +
    labs(y = "View", x = element_blank(), subtitle = 'Factor 3') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)), axis.text.y = element_text(size = rel(1)))
)
# 3.78x3.16

wrap_plots(pH, (pI2 + pI1 + plot_layout(guides = "collect"))) + plot_layout(design = "1222222")

