library(magrittr)
library(dplyr)
exposure <- read.csv("data/AirExposure_CensusTract_OCEAN_Mapping.csv",
                 row.names = 1) %>%
  tibble::rownames_to_column("EdgarID") %>%
  dplyr::mutate(across(everything(), ~ dplyr::na_if(., "."))) %>%
  dplyr::mutate(pm_2_5 = pm_2.5) %>%
  dplyr::select(EdgarID, pm_2_5, bc, so4, starts_with("Contact")) %>%
  dplyr::mutate(across(c(pm_2_5, bc, so4), as.numeric))

id_map <- readxl::read_xlsx("data/MapFile_For_Charlotte.xlsx")


clinical <- readxl::read_xlsx("data/Heidelberg_Neptune_clinical_N120_20240513.xlsx") %>%
  dplyr::left_join(id_map)

modules <- readxl::read_xlsx("data/OCEAN_NEPTUNE_SAFIDs_Modules.xlsx") %>%
  dplyr::left_join(id_map)

metadatav7 <- read.csv("data/OCEAN_Metadata_v7.csv") %>%
  dplyr::mutate(across(everything(), ~na_if(., "."))) %>% 
  dplyr::mutate(pm_2_5 = pm_2.5) %>%
  dplyr::select(-pm_2.5) %>%
  dplyr::rename_with(~ gsub("MEST.C", "MESTC", .x, fixed = TRUE)) %>%
  dplyr::rename_with(~ gsub(".", "_", .x, fixed = TRUE)) %>%
  dplyr::mutate(across(c(Age, eGFR_at_Bx_NEPTUNE, UPCR_at_Bx_NEPTUNE,
                         number_of_APOL1_risk_alleles_WGS1_2, 
                         eGFRatBL_NEPTUNE, UPCRatBL_NEPTUNE,
                         TNF_Grouping_JASN, nGloms_NEPTUNE, SMHpct_NEPTUNE,
                         GMHpct_NEPTUNE, SEHpct_NEPTUNE, GEHpct_NEPTUNE,
                         MESHYPERpct_NEPTUNE, SEGEPILESpct_NEPTUNE,
                         GLOEPILESpct_NEPTUNE, CELLCRSpct_NEPTUNE, 
                         FIBROCELLCRSpct_NEPTUNE,
                         InterstitialFibrosis_NEPTUNE, TubularAtrophy_NEPTUNE,
                         pm_2_5, bc, so4,
                         MESTC_nGloms, MESTC_nGSG, 
                         MESTC_nGlom_c_MC_hypercell,
                         MESTC_nGlom_c_Endo_hypercell, MESTC_nSSG,
                         MESTC_perCent_TA_IF_Cortex, 
                         MESTC_nGlom_Cell_Fibrocell_Crescent,
                         MESTC_TotalScore), ~as.numeric(.))) %>%
  dplyr::mutate(log_UPCR_at_Bx_NEPTUNE = log(UPCR_at_Bx_NEPTUNE))

# metadata <- read.csv("results/preprocessing/metadata_clean.csv", row.names = 1) %>%
#   dplyr::select(-c(pm_2_5, bc, so4, starts_with("Contact"))) %>%
#   tibble::rownames_to_column("EdgarID")

updated_metadata <- metadatav7 %>%
  dplyr::left_join(clinical) %>%
  dplyr::left_join(modules) %>%
  dplyr::left_join(exposure) %>%
  dplyr::distinct() %>%
  dplyr::select_if(~!all(is.na(.))) %>%
  dplyr::mutate(across(everything(), ~stringr::str_replace(., "[:digit:]: ", "")))

# metadata <- read.csv("results/preprocessing/metadata_clean.csv",
#                      row.names = 1) %>%
#   tibble::rownames_to_column("EdgarID")

# metadata <- metadata %>%
#   dplyr::full_join(data)

# write.csv(metadata, "results/preprocessing/metadata_clean.csv",
#           row.names = FALSE)

write.csv(updated_metadata, "results/preprocessing/updated_metadata_clean.csv",
          row.names = FALSE)

test <- read.csv("results/preprocessing/metadata_clean.csv",
                 row.names = 1) %>%
  tibble::rownames_to_column("EdgarID")
