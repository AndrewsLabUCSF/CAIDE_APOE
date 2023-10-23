setwd("/Users/elinorvelasquez/Desktop/adni_nacc/")
library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("readxl")
library("readr")
library("webshot2")
library("ggplot2")
library("gtsummary")
library("broom")
library("MASS")

########################################################
# Read data from the combined adni-nacc datasets
########################################################
# Case 1:
data_adni_nacc <- read_csv("~/Desktop/adni_nacc/wynton/adni_nacc.csv") %>%
  janitor::clean_names()

# was: caide_1
# Shea: tbl1_df <- data_wo_chol_case1 %>%
caide_case1 <- data_adni_nacc %>%
  mutate(
    Gender = gender,
    diag = case_when(
        dx_bl == 'AD' ~ "Case",
        dx_bl == 'ADRD' ~ "Case",
        dx_bl == 'CN' ~ "Control",
        dx_bl == 'EMCI' ~ "Control",
        dx_bl == 'LMCI' ~ "Case",
        dx_bl == "MCI" ~ "Case",
        dx_bl == 'SMC' ~ "Control"
    ),
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
  ) 

adni_nacc_case1 <- caide_case1  %>%
  mutate(
    caide_apoe1 = fct_relevel(caide_apoe1, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe1 = fct_relevel(mcaide_apoe1, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+")
  ) %>%
  dplyr::select(race, diag,
                #Hypercholesterolemia3, 
                cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                caide1, apoe, caide_apoe1, mcaide1, mcaide_apoe1 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", diag ~ "categorical",
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
                 race ~ "Race", diag ~ "Diagnosis", cohort ~ "Cohort",
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
    filename = "~/Desktop/adni_nacc/table_adni_nacc_1.png", 
    vwidth = 648, vheight = 288
  )

# webshot("~/Desktop/adni_nacc/table_adni_nacc_wo_chol_case1.html", 
#        "~/Desktop/adni_nacc/table_adni_nacc_wo_chol_case1.png")

############################################################
# Case 3. m/caide table w/ imputed chol
############################################################

adni_nacc_case3 <- caide_case1 %>%
  mutate(
    Hypercholesterolemia3 = case_when(
      (hypchol3 == 1) ~ "Yes",
      (hypchol3 == 0) ~ "No",
      is.na(hypchol3) ~ NA
    ),
    caide_apoe3 = fct_relevel(caide_apoe3, "Low, e2+", "Low, e3/e3", "Low, e4+",
                                "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                                "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe3 = fct_relevel(mcaide_apoe3, "Low, e2+", "Low, e3/e3", "Low, e4+",
                                 "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                                 "High, e2+", "High, e3/e3", "High, e4+")
    ) %>%
  dplyr::select(race, diag, cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                Hypercholesterolemia3,
                caide3, apoe, caide_apoe3, mcaide3, mcaide_apoe3 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", diag ~ "categorical",
                Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", race ~ "categorical",
                caide3 ~ "continuous", caide_apoe3 ~ "categorical",
                mcaide3 ~ "continuous", mcaide_apoe3 ~ "categorical",
                Hypercholesterolemia3 ~ "categorical"
    ),  
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide3 ~ "CAIDE_imputed_chol", mcaide3 ~ "mCAIDE_imputed_chol",
                 apoe ~ "APOE", caide_apoe3 ~ "CAIDE x APOE_imputed_chol",
                 mcaide_apoe3 ~ "mCAIDE x APOE_imputed_chol",
                 Hypercholesterolemia3 ~ "Hypercholesterolemia"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave("table_adni_nacc_3.png", vwidth = 576, vheight = 288)

############################################################
# Case 2. m/caide table w/ chol: data_w_chol_noNA_case2
############################################################

adni_nacc_case2 <- caide_case1 %>%
  mutate(
    Hypercholesterolemia2 = case_when(
      (hypchol2 == 1) ~ "Yes",
      (hypchol2 == 0) ~ "No",
      is.na(hypchol2) ~ NA
    ),
    caide_apoe2 = fct_relevel(caide_apoe2, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe2 = fct_relevel(mcaide_apoe2, "Low, e2+", "Low, e3/e3", "Low, e4+",
                               "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                               "High, e2+", "High, e3/e3", "High, e4+")
  ) %>%
  dplyr::select(race, diag, cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                Hypercholesterolemia2,
                caide2, apoe, caide_apoe2, mcaide2, mcaide_apoe2 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", diag ~ "categorical",
                Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", race ~ "categorical",
                caide2 ~ "continuous", caide_apoe2 ~ "categorical",
                mcaide2 ~ "continuous", mcaide_apoe2 ~ "categorical",
                Hypercholesterolemia2 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide2 ~ "CAIDE_w_cholNAs", mcaide2 ~ "mCAIDE_w_cholNAs",
                 apoe ~ "APOE", caide_apoe2 ~ "CAIDE x APOE_w_cholNAs",
                 mcaide_apoe2 ~ "mCAIDE x APOE_w_cholNAs",
                 Hypercholesterolemia2 ~ "Hypercholesterolemia"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave("table_adni_nacc_2.png", vwidth = 576, vheight = 288)

########################################################
# 2. mcaide table w/o chol --OLD CODE--
########################################################

mcaide_1 <- data_wo_chol %>%
  mutate(
    Gender = case_when(
      (caide_sex == 1) ~ "Male",
      (caide_sex == 0) ~ "Female",
      is.na(caide_sex) ~ NA
    ),
    Age = case_when(
      (i_age < 64) ~ "< 64",
      (i_age >= 64 & i_age <= 73) ~ "64 - 73",
      (i_age > 73) ~ "> 73",
      is.na(i_age) ~ NA
    ),
    Education = case_when(
      (i_education < 12) ~ "< 12",
      (i_education >= 12 & i_education <= 16) ~ "12 - 16",
      (i_education > 16) ~ "> 16", is.na(i_education) ~ NA
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
    )
  )

mcaide_1 %>%
  dplyr::select(
    race,
    Gender, Age, Education, Obesity, Hypertension, 
    mcaide, apoe, mcaide_apoe, 
    dx_bl) %>%
  tbl_summary(
    by = race,
    type = list(mcaide ~ "continuous", Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", 
                race ~ "categorical",
                dx_bl ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list( 
      dx_bl ~ "Diagnosis", Gender ~ "Gender", Age ~ "Age", 
      Education ~ "Education", 
      Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
      mcaide ~ "MCAIDE", race ~ "Race",
      apoe ~ "APOE", mcaide_apoe ~ "MCAIDE x APOE"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "table_mcaide_wo_chol_test.png", vwidth = 648, vheight = 288
  )

# webshot("table_mcaide_wo_chol.html", "table_mcaide_wo_chol.png")

############################################################
# 1. caide table w/ chol
############################################################

caide_1 <- data_w_chol %>%
  mutate(
    Gender = case_when(
      (caide_sex == 1) ~ "Male",
      (caide_sex == 0) ~ "Female",
      is.na(caide_sex) ~ NA
    ),
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
    ),
    Hypercholesterolemia = case_when(
      (hypchol == 1) ~ "Yes",
      (hypchol == 0) ~ "No",
      is.na(hypchol) ~ NA
    )
  )

caide_1 %>%
  dplyr::select(dx_bl, race,
                Gender, Age, Education, Obesity, Hypertension, 
                caide, apoe, caide_apoe, Hypercholesterolemia 
  ) %>%
  tbl_summary(
    by = race,
    type = list(caide ~ "continuous", Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", race ~ "categorical",
                dx_bl ~ "categorical", Hypercholesterolemia ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", dx_bl ~ "Diagnosis",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide ~ "CAIDE", 
                 apoe ~ "APOE", caide_apoe ~ "CAIDE x APOE",
                 Hypercholesterolemia ~ "Hypercholesterolemia"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "table_caide_w_chol_test.png", vwidth = 648, vheight = 288
  )

# webshot("table_caide_w_chol.html", "table_caide_w_chol.png")

########################################################
# 2. mcaide table w/ chol
########################################################

mcaide_1 <- data_w_chol %>%
  mutate(
    Gender = case_when(
      (caide_sex == 1) ~ "Male",
      (caide_sex == 0) ~ "Female",
      is.na(caide_sex) ~ NA
    ),
    Age = case_when(
      (i_age < 64) ~ "< 64",
      (i_age >= 64 & i_age <= 73) ~ "64 - 73",
      (i_age > 73) ~ "> 73",
      is.na(i_age) ~ NA
    ),
    Education = case_when(
      (i_education < 12) ~ "< 12",
      (i_education >= 12 & i_education <= 16) ~ "12 - 16",
      (i_education > 16) ~ "> 16", is.na(i_education) ~ NA
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
    ),
    Hypercholesterolemia = case_when(
      (hypchol == 1) ~ "Yes",
      (hypchol == 0) ~ "No",
      is.na(hypchol) ~ NA
    )
  )

mcaide_1 %>%
  dplyr::select(
    race,
    Gender, Age, Education, Obesity, Hypertension, 
    mcaide, apoe, mcaide_apoe,
    Hypercholesterolemia, dx_bl) %>%
  tbl_summary(
    by = race,
    type = list(mcaide ~ "continuous", Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", 
                race ~ "categorical",
                dx_bl ~ "categorical",
                Hypercholesterolemia ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list( 
      dx_bl ~ "Diagnosis", Gender ~ "Gender", Age ~ "Age", 
      Education ~ "Education", 
      Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
      mcaide ~ "MCAIDE", race ~ "Race",
      apoe ~ "APOE", mcaide_apoe ~ "MCAIDE x APOE",
      Hypercholesterolemia ~ "Hypercholesterolemia"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "table_mcaide_w_chol_test.png", vwidth = 648, vheight = 288
    ) # 9 x 4 inches

# webshot("table_mcaide_w_chol.html", "table_mcaide_w_chol.png")
