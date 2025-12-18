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

my_comparisons <- list(c("None", "Pre V3"), c("None", "At V3"))
ps1 <- df %>% 
  dplyr::mutate(Steroid = factor(ifelse(
    SteroidAtV3_NEPTUNE == "Yes",
    "At V3",
    ifelse(PAT_SteroidPreV3_NEPTUNE == "Yes",
           "Pre V3", ifelse(
             PAT_SteroidPreV3_NEPTUNE == "No",
             "None", NA
           ))),levels = c("None", "Pre V3", "At V3"))) %>%
  tidyr::drop_na(Steroid) %>%
  ggplot(aes(x = Steroid, y = Factor5)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Age), size= 1) +
  scale_color_viridis_c(direction = -1,limits = c(0,80)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(subtitle = "Adult and pediatric\nparticipants") + 
  ylim(-0.2, 0.8) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons

ps2 <- df %>% 
  dplyr::filter(Pediatric_NEPTUNE == "Yes") %>%
  dplyr::mutate(Steroid = factor(ifelse(
    SteroidAtV3_NEPTUNE == "Yes",
    "At V3",
    ifelse(PAT_SteroidPreV3_NEPTUNE == "Yes",
           "Pre V3", ifelse(
             PAT_SteroidPreV3_NEPTUNE == "No",
             "None", NA
           ))),levels = c("None", "Pre V3", "At V3"))) %>%
  tidyr::drop_na(Steroid) %>%
  ggplot(aes(x = Steroid, y = Factor5)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Age), size= 1) +
  scale_color_viridis_c(direction = -1, limits = c(0,80)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(subtitle = "Pediatric participants") +
  ylim(-0.2, 0.8) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons

ps3 <- df %>%
  dplyr::filter(Pediatric_NEPTUNE == "No") %>%
  dplyr::mutate(Steroid = factor(
    ifelse(
      SteroidAtV3_NEPTUNE == "Yes",
      "At V3",
      ifelse(
        PAT_SteroidPreV3_NEPTUNE == "Yes",
        "Pre V3",
        ifelse(PAT_SteroidPreV3_NEPTUNE == "No", "None", NA)
      )
    ),
    levels = c("None", "Pre V3", "At V3")
  )) %>%
  tidyr::drop_na(Steroid) %>%
  ggplot(aes(x = Steroid, y = Factor5)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = Age), size = 1) +
  scale_color_viridis_c(direction = -1, limits = c(0, 80)) +     theme_classic() +  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  labs(subtitle = "Adult participants") + ylim(-0.2, 0.8) +
  ggpubr::stat_compare_means(comparisons = list(c("None", "At V3")), method = "t.test") # Add pairwise comparisons

ps1 + ps2 + ps3 + patchwork::plot_layout(guides = "collect")
