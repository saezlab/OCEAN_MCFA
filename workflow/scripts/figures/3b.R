# Figure 3b: Number of patients
library(ggplot2)
library(scales)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(magrittr)
library(patchwork)
library(ggpubr)

metadata <- read.csv("data/metadata.csv")
factors <- read.csv('results/mofacell/widescores~combined_sn10xLRCluster2.csv')

df <- metadata %>%
  dplyr::filter(Prep == "SN" & Method == "10x") %>%
  #dplyr::filter(OCEAN_MS_Exp_ID %in% factors$OCEAN_MS_Exp_ID) %>%
  dplyr::mutate(Cohort_Manuscript = factor(Cohort_Manuscript,
                                           levels = c("FSGS", "IgAN", "MCD", "OtherProteinuric", "NotApplicable", "LD", "Transplant"),
                                           labels = c("FSGS", "IgAN", "MCD", "Other Proteinuric", "Not Applicable", "Living Donor", "Transplant"))) %>%
  dplyr::filter(Cohort_Manuscript != "Not Applicable") %>%
  dplyr::select(EdgarID, Cohort, Cohort_Manuscript) %>%
  dplyr::count(Cohort_Manuscript)

(
  p <- ggplot(df, aes(x = Cohort_Manuscript, y = n)) +
    geom_col() +
    theme_classic() + 
    labs(subtitle = "10x snRNA-seq", x = element_blank(), y = "Number of samples") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
)
