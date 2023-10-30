setwd("/Users/elinorvelasquez/Desktop/ADNI/")

library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("broom")

adni_raw <- read_csv("~/Desktop/ADNI/DXSUM_PDXCONV_ADNIALL_12Oct2023.csv") %>%
  janitor::clean_names()

adni_novel <- adni_raw %>%
  mutate(
    diag = case_when(
      ### ADRD ###
      dxcurren == 3 | dxad == 1 | dxddue == 1 ~ 1, # all cohorts; AD
      
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") & 
        dxddue == 2 & dxodes == 1 ~ 1, # Frontal-temporal
      phase == "ADNI1" & dxothdem == 1 & dxodes == 1 ~ 1, # Frontal-temporal
      
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") & 
        dxddue == 2 & dxodes == 6 ~ 1, # NPH
      phase == "ADNI1" & dxothdem == 1 & dxodes == 4 ~ 1, # NPH
      
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") & 
        dxddue == 2 & dxodes == 9 ~ 1, # Vascular Dementia
      phase == "ADNI1" & dxothdem == 1 & dxodes == 7 ~ 1, # Vascular Dementia
      
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") & 
        dxddue == 2 & dxodes == 8 ~ 1, # Corticobasal Degeneration
      
      (phase == "ADNIGO" | phase == "ADNI2") & 
        dxddue == 2 & dxodes == 13 ~ 1, # Posterior Cortical Dysfunction
      phase == "ADNI3" & dxddue == 2 & dxodes == 13 ~ 1, # Posterior Cortical Dysfunction
      phase == "ADNI1" & dxothdem == 1 & dxodes == 11 ~ 1, # Posterior Cortical Dysfunction
      
      ### CN ###
      dxcurren == 1 | # ADNI1
      dxnorm == 1 | # ADNI1
      diagnosis == 1 | # ADNI3
      dxchange == 1 ~ 0, # ADNIGO, ADNI2
      
      ### MCI ###
      diagnosis == 2 | # ADNI3
      dxcurren == 2 | # ADNI1
      dxmci == 1 | # ADNI1
      dxchange == 2 | # ADNIGO, ADNI2
      dxmdue == 1 ~ 1, # ADNIGO, ADNI2, ADNI3
      phase == "ADNI1" & 
        (dxmothet == 1 | dxmothet == 5 | dxmothet == 6) ~ 1,
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") &
        (dxmothet == 1 | dxmothet == 8 | dxmothet == 9) ~ 1,
      dxmothsp == "Dementia with Lewy bodies" |
        dxmothsp == "vascular cognitive impairment" |
        dxmothsp == "'Mixed' Vascular/Alzheimer's" |
        dxmothsp == "(Also Alzheimer's Disease)" |
        dxmothsp == "MCI due to Alzheimer's Disease" |
        dxmothsp == "Mixed: Vascular/Alzheimer's" | 
        dxmothsp == "vascular dis" |
        dxmothsp == "Alzheimer's Disease" |
        dxmothsp == "Dementia with Lewy Bodies" | 
        dxmothsp == "Lewy Bodies" |
        dxmothsp == "Lewy Body Disease" |
        dxmothsp == "Mixed AD and vascular cognitive impairment. The amount of 
      vascular disease cannot be ignored in formulating the dianosis" |
        dxmothsp == "Possible Limbic associate TAU Encephalopathy (LATE)" |
        dxmothsp == "TDP-43" |
        dxmothsp == "Tauopathy" | 
        dxmothsp == "also due to Alzheimer's disease" |
        dxmothsp == "fibromyalgia (box is also checked for MCI due to AD at the 
      top on source doc)" |
        dxmothsp == "hydrocephalus" | # NPH
        dxmothsp == "mixed AD and Vascular" |
        dxmothsp == "mixed AD and vascular" |
        dxmothsp == "mixed Alz and vascular" |
        dxmothsp == "mixed etiology = MCI due to AD and MCI due to Vascular 
      Dementia" | 
        dxmothsp == "suspected mixed etiology (AD, vascular, LATE)" |
        dxmothsp == "vascular MCI" |
        dxmothsp == "Possible Lewy Body Disease vs. Alzheimer Disease" |
        dxmothsp == "aging or early AD - cannot tell for sure" ~ 1,
      
        dxmothsp == "-4" |
        dxad == "-4" |
        dxodes == "-4" |
        dxothdem == "-4" |
        dxnorm == "-4" |
        dxmci == "-4" |
        dxmdue == "-4" |
        dxmothet == "-4" ~ NA_real_,
        ###TRUE ~
    )
  )

