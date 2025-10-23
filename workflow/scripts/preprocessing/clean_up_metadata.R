# TODO: write rule
library(magrittr)

print("INFO: Job running...")
if (exists('snakemake')){
  metadata <- read.csv(snakemake@input$metadata, row.names = 1)
  metadata_clean_csv <- snakemake@output$metadata_clean
} else {
  metadata <- read.csv('data/Metadata_For_Charlotte.csv', row.names = 1)
  metadata_clean_csv <- 'results/preprocessing/metadata_clean.csv'
}

print(head(metadata))
simplify_cohort <- function(cohort){
  ifelse(stringr::str_detect(cohort, "IgAN"), "IgAN",
         ifelse(cohort %in% c("AIN/T2D", "C1QN", "C3_Glomerulopathy",
                              "COVAN", "Immune_Cmpx_GN",
                              "MN_Vax", "MPGN"), "Other",
                ifelse(stringr::str_detect(cohort, "FSGS"), "FSGS", cohort)))
}

metadata_clean <- metadata %>%
  dplyr::mutate(across(where(is.character), ~dplyr::na_if(., "N/A"))) %>%
  dplyr::mutate(across(c(Age, eGFRatBx_NEPTUNE, UPCRatBx_NEPTUNE,
                         number_of_APOL1_risk_alleles, nGloms, SMHpct,
                         GMHpct, SEHpct, GEHpct, MESHYPERpct, SEGEPILESpct,
                         GLOEPILESpct, CELLCRSpct, FIBROCELLCRSpct,
                         InterstitialFibrosis, TubularAtrophy), ~as.numeric(.))) %>%
  dplyr::mutate(Cohort_Charlotte = simplify_cohort(Cohort))

write.csv(metadata_clean, metadata_clean_csv)
