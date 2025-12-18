library(magrittr)
library(ggplot2)
library(ggpubr)

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
                "log(uPCR)" = log(NEPTUNE_UPCR_at_Bx))

df %>% 
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
  dplyr::mutate(eGFRAgeStratification = ifelse(Age < 18 & eGFR >= 45,
                                               "YoungHigher",
                                               ifelse(Age >= 18 & eGFR < 45,
                                                      "AdultLow", "Other"))) %>%
  dplyr::count(Group, eGFRAgeStratification) %>%
  dplyr::group_by(eGFRAgeStratification) %>% dplyr::mutate(totaleGFRAge = sum(n))

pF <- df %>% 
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
  dplyr::mutate(Group = factor(Group, levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5",   "HighF3HighF5"))) %>%
  ggplot(aes(x = Age, y = eGFR)) + geom_hline(yintercept = 45, alpha = 0.3) + geom_vline(xintercept = 18, alpha = 0.3) + geom_point(aes(color = Group)) + theme_classic() + scale_color_brewer(palette="Set1") + ylim(0, 180) +
  labs(y = "eGFR at biopsy (ml/min/1.73 m²)", x = "Age (y)")

pF

my_comparisons <- list(c("LowF3HighF5", "LowF3LowF5"),
                       c("LowF3LowF5", "HighF3LowF5"),
                       c("LowF3HighF5", "HighF3LowF5"))
(
  pG <- df %>%
    dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                                 ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                        ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                               ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
    dplyr::mutate(Group = factor(Group, levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5",   "HighF3HighF5"))) %>%
    dplyr::filter(Group %in% c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5")) %>% 
    ggplot(aes(x = forcats::fct_rev(Group),
               y = eGFR)) + 
    geom_boxplot() + geom_jitter(aes(color = Group), size = 1) + theme_classic() + scale_color_brewer(palette="Set1") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + # Add pairwise comparisons
    labs(x = element_blank()) + ylim(0,220) +
    ggpubr::stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") + labs(y = "eGFR at biopsy (ml/min/1.73 m²)")
)

# Supplementary Figure
my_comparisons <- list(c("LowF3HighF5", "LowF3LowF5"),
                       c("LowF3LowF5", "HighF3LowF5"),
                       c("LowF3HighF5", "HighF3LowF5"))
p8 <- df %>%
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
  dplyr::mutate(Group = factor(Group, levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5",   "HighF3HighF5"))) %>%
  dplyr::filter(Age >= 18 & eGFR >= 45 & Group %in% c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5")) %>% 
  ggplot(aes(x = forcats::fct_rev(Group),
             y = eGFR)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = Group), height = 0) + 
  theme_classic() + 
  scale_color_brewer(palette="Set1") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test") + # Add pairwise comparisons 
  labs(y = "eGFR at biopsy (ml/min/1.73 m²)",
       x = element_blank(),
       subtitle = "Age ≥ 18,\neGFR ≥ 45 ml/min/1.73 m²")  +
  ylim(0,200) 

p9 <- df %>%
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
  dplyr::mutate(Group = factor(Group, levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5",   "HighF3HighF5"))) %>%
  dplyr::filter(Age >= 18 & Group %in% c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5")) %>% 
  ggplot(aes(x = forcats::fct_rev(Group),
             y = eGFR)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = Group), height = 0) + 
  theme_classic() + 
  scale_color_brewer(palette="Set1") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test") + # Add pairwise comparisons
  labs(y = "eGFR at biopsy (ml/min/1.73 m²)",
       x = element_blank(),
       subtitle = "Age ≥ 18")  +
  ylim(0,200) 

p10 <- df %>%
  dplyr::mutate(Group = ifelse(Factor3 < 0 & Factor5 > 0, "LowF3HighF5",
                               ifelse(Factor3 > 0 & Factor5 < 0, "HighF3LowF5",
                                      ifelse(Factor3 < 0 & Factor5 < 0, "LowF3LowF5", 
                                             ifelse(Factor3 > 0 & Factor5 > 0, "HighF3HighF5", NA))))) %>%
  dplyr::mutate(Group = factor(Group, levels = c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5",   "HighF3HighF5"))) %>%
  dplyr::filter(Age < 18 & Group %in% c("HighF3LowF5", "LowF3LowF5", "LowF3HighF5")) %>% 
  ggplot(aes(x = forcats::fct_rev(Group),
             y = eGFR)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = Group), height = 0) + 
  theme_classic() + 
  scale_color_brewer(palette="Set1") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test") + # Add pairwise comparisons 
  labs(y = "eGFR at biopsy (ml/min/1.73 m²)",
       x = element_blank(),
       subtitle = "Age < 18") +
  ylim(0,200) 

p9 + p10
