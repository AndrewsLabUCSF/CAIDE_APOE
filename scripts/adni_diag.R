setwd("/Users/elinorvelasquez/Desktop/ADNI/")

library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("broom")
library("webshot2")
library("reticulate")

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
        dxddue == 2 & dxodes == 9 ~ 1, # Vascular Dementia
      phase == "ADNI1" & dxothdem == 1 & dxodes == 7 ~ 1, # Vascular Dementia
      
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
        (dxmothet == 1 | # Frontal Lobe Dementia
           dxmothet == 6) ~ 1, # Vascular Dementia
      (phase == "ADNIGO" | phase == "ADNI2" | phase == "ADNI3") &
        (dxmothet == 1 | # Frontal Lobe Dementia
           dxmothet == 9) ~ 1, # Vascular Dementia
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
        TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::filter(!is.na(diag))    

######################################
# Make table with this new diagnosis variable
######################################

adni_nacc_raw <- read_csv("~/Desktop/adni_nacc/wynton/adni_nacc.csv") %>%
  janitor::clean_names()

janitor::get_dupes(adni_nacc_raw) # no dupes

adni_nacc_raw %>% count(cohort)

# Left-join new diagnosis to adni_nacc data:

adni_novel$rid = as.factor(adni_novel$rid)

adni_nacc_diag <- adni_nacc_raw %>%
  mutate(
    diag = case_when(
      dx_bl == "AD" ~ 1,
      dx_bl == "LMCI" ~ 1,
      dx_bl == "EMCI" ~ 0,
      dx_bl == "CN" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::filter(!is.na(diag))

adni_nacc_diag %>% count(diag) 

adni_nacc_novel <- adni_nacc_diag %>% 
  left_join(adni_novel, by=c('rid', 'diag')) %>%
  dplyr::filter(!is.na(diag))

adni_nacc_novel %>% count(cohort)
adni_nacc_novel %>% count(diag)

adni_nacc_novel %>% janitor::get_dupes() # no dupes

# Make new variables for the table: Case 1:

table_case1 <- adni_nacc_novel %>%
  mutate(
    Gender = gender,
    Diag = case_when(
      diag == 1 ~ "ADRD, MCI",
      diag == 0 ~ "Cog Normal",
      TRUE ~ NA_character_
    ),
    #diag = case_when(
    #  dx_bl == 'AD' ~ "Case",
    #  dx_bl == 'ADRD' ~ "Case",
    #  dx_bl == 'CN' ~ "Control",
    #  dx_bl == 'EMCI' ~ "Control",
    #  dx_bl == 'LMCI' ~ "Case",
    #  dx_bl == "MCI" ~ "Case",
    #  dx_bl == 'SMC' ~ "Control"
    #),
    Age = case_when(
      (i_age < 47) ~ "< 47",
      (i_age >= 47 & i_age <= 53) ~ "47 - 53",
      (i_age > 53) ~ "> 53",
      is.na(i_age) ~ NA
    ),
    Education = case_when(
      (i_education < 7) ~ "< 7",
      (i_education >= 7 & i_education <= 9) ~ "7 - 9",
      (i_education > 9) ~ "> 9", is.na(i_education) ~ NA
    ),
    Obesity = case_when(
      (i_bmi <= 30) ~ "No",
      (i_bmi > 30) ~ "Yes",
      is.na(i_bmi) ~ NA
    ),
    Hypertension = case_when(
      (i_systbp <= 140) ~ "No",
      (i_systbp > 140) ~ "Yes",
      is.na(i_systbp) ~ NA
    )#,
    #Hypercholesterolemia3 = case_when(
    #  (hypchol3 == 1) ~ "Yes",
    #  (hypchol3 == 0) ~ "No",
    #  is.na(hypchol3) ~ NA
    #)
  ) # end of caide_case1

table_case1 %>% count(Diag)

adni_nacc_case1 <- table_case1  %>%
  mutate(
    caide_apoe1 = fct_relevel(caide_apoe1, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe1 = fct_relevel(mcaide_apoe1, "Low, e2+", "Low, e3/e3", "Low, e4+",
                               "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                               "High, e2+", "High, e3/e3", "High, e4+")
  ) %>%
  dplyr::select(race, Diag,
                #Hypercholesterolemia3, 
                cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                caide1, apoe, caide_apoe1, mcaide1, mcaide_apoe1 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", Diag ~ "categorical",
                Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", race ~ "categorical",
                caide1 ~ "continuous", caide_apoe1 ~ "categorical",
                mcaide1 ~ "continuous", mcaide_apoe1 ~ "categorical"#,
                #Hypercholesterolemia3 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", Diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide1 ~ "CAIDE_no_chol", mcaide1 ~ "mCAIDE_no_chol",
                 apoe ~ "APOE", caide_apoe1 ~ "CAIDE x APOE_no_chol",
                 mcaide_apoe1 ~ "mCAIDE x APOE_no_chol"#,
                 #Hypercholesterolemia3 ~ "Hypercholesterolemia_imp_chol"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "~/Desktop/adni_nacc/table_adni_nacc_1.html", 
  )

webshot("~/Desktop/adni_nacc/table_adni_nacc_1.html", 
        "~/Desktop/adni_nacc/table_adni_nacc_1.png",
        vwidth = 648, vheight = 288)










