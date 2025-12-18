library(magrittr)
library(ggplot2)
library(ggpubr)

pod_psbulk <- read.csv("results/mofacell/pod_psbulk.csv", row.names = 1) %>%
  tibble::rownames_to_column('OCEAN_MS_Exp_ID')

metadata <- read.csv("data/metadata.csv")
scores <- read.csv("results/mofacell/widescores~combined_sn10xLRCluster2.csv") %>%
  dplyr::mutate(Factor5=-Factor5)

df <- dplyr::left_join(scores, metadata, by = c("OCEAN_MS_Exp_ID")) %>%
  dplyr::mutate("Project (Batch)" = Project1,
                Age,
                "Cohort" = Cohort_Manuscript,
                "Self-ID Race" = Race,
                "eGFR" = NEPTUNE_eGFR_at_Bx,
                "Steroid at V3" = SteroidAtV3_NEPTUNE,
                "Interstitial Fibrosis" = InterstitialFibrosis_NEPTUNE,
                "Tubular Atrophy" = TubularAtrophy_NEPTUNE,
                "log(uPCR)" = log(NEPTUNE_UPCR_at_Bx)) %>%
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) 
# 
apol1_POD <- pod_psbulk %>%
  dplyr::transmute(OCEAN_MS_Exp_ID, APOL1) %>%
  dplyr::left_join(df, by = 'OCEAN_MS_Exp_ID') %>%
  dplyr::mutate(Group = factor(
    Group,
    levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5", "HighF3HighF5")
  ))

my_comparisons <- list(c("HighF3LowF5", "LowF3HighF5"))


p3 <- apol1_POD %>%
  dplyr::filter(Cohort %in% c("FSGS")) %>%
  dplyr::filter(Group %in% c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5")) %>% 
  ggplot(aes(x = forcats::fct_rev(Group),
             y = APOL1)) + 
  geom_boxplot() + geom_jitter(aes(color = Group)) + theme_classic() +
  scale_color_brewer(palette="Set1") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + # Add pairwise comparisons
  labs(x = element_blank()) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", method.args = list(alternative = "greater")) +
  labs(subtitle = "FSGS:POD", y = "log1p normalised APOL1 expr.")  + ylim(c(0,1.5))

p3
