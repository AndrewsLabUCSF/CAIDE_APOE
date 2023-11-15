
library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("broom")
library("ggplot2")
library("gtsummary")
library("rmarkdown")
library(vcfR)

setwd("~/gitcode/CAIDE_APOE")

## Import datasets
### Key variables 
merge.raw <- read_csv("resources/ADNI/ADNIMERGE_25Sep2023.csv") %>%
  janitor::clean_names()

### Diagnostic summary
dxsum_raw <- read_csv("resources/ADNI/DXSUM_PDXCONV_ADNIALL_25Sep2023.csv", guess_max = 14000) %>%
  janitor::clean_names() 

### Vitals 
vitals.raw <- read_csv("resources/ADNI/VITALS_26Sep2023.csv") %>%
  janitor::clean_names()

### Cholesterol 
labs.raw <- read_csv("resources/ADNI/LABDATA_26Sep2023.csv") %>%
  janitor::clean_names() 

### Genetics
apoeres.raw <- read_csv("resources/ADNI/APOERES_15Nov2023.csv") %>%
  janitor::clean_names() 

## Wrangling 
### Genetics
apoeres <- apoeres.raw %>%
  unite(apoe_geno, apgen1, apgen2, sep = "") %>%
  select(phase, rid, ptid, viscode, apoe_geno)

### Merge 
merge <- merge.raw %>%
  mutate(
    race = case_when(
      # ptraccat != "Asian" & ptethcat != "Hisp/Latino" & ptraccat != "Black" 
      #  & ptraccat != "White" ~ "Other",
      ptraccat == "White" & ptethcat != "Hisp/Latino" ~ "Non-Hispanic White",
      ptraccat == "White" & ptethcat == "Hisp/Latino" ~ "Hispanic",
      ptraccat == "More than one" & ptethcat == "Hisp/Latino" ~ "Hispanic",
      ptraccat == "Black" ~ "Black",
      ptraccat == "Asian" ~ "Asian",
      ptraccat == "Unknown" ~ "Other",
      ptraccat == "Am Indian/Alaskan" ~ "Other",
      ptraccat == "Hawaiian/Other PI" ~ "Other",
      ptraccat == "More than one" ~ "Other",
      TRUE ~ NA_character_
    )) %>%
  filter(viscode == "bl") %>%
  select(rid, ptid, origprot, colprot, viscode, examdate_bl, examdate, cdrsb, 
         ptraccat, ptethcat, age, ptgender, pteducat, apoe4, race)

### ADRD Diagnosis aligned with NACC UDS
dxsum <- dxsum_raw %>%
  mutate(
    naccudsd = case_when(
      dxnorm == 1 ~ 1, 
      dxmci == 1 ~ 3, 
      dxad == 1 ~ 4, 
      dxothdem == 1 ~ 4, 
      dxchange == 1 ~ 1, 
      dxchange == 2 ~ 3, 
      dxchange == 3 ~ 4, 
      dxchange == 4 ~ 3, 
      dxchange == 5 ~ 4, 
      dxchange == 6 ~ 4, 
      dxchange == 7 ~ 1, 
      dxchange == 8 ~ 3, 
      dxchange == 9 ~ 1, 
      diagnosis == 1 ~ 1, 
      diagnosis == 2 ~ 3, 
      diagnosis == 3 ~ 4, 
    ), 
    naccetpr = case_when(
      naccudsd == 1 ~ 88,
      dxmdue == 1 ~ 1,
      dxddue == 1 ~ 1,
      dxad == 1 ~ 4, 
      phase == "ADNI1" & dxmdue == 1 ~ 1, 
      phase == "ADNI1" & dxmothet == 1 ~ 7, 
      phase == "ADNI1" & dxmothet == 6 ~ 8, 
      phase == "ADNI1" & dxmothet == "6:08" ~ 8, 
      phase == "ADNI1" & dxmothet == 8 ~ 30, 
      phase == "ADNI1" & dxodes == 1 ~ 7, 
      phase == "ADNI1" & dxodes == 12 ~ 18, 
      dxmothet == 1 ~ 7,
      dxmothet == 2 ~ 18,
      dxmothet == 3 ~ 11,
      dxmothet == 4 ~ 4,
      dxmothet == 5 ~ 26,
      dxmothet == 6 ~ 14,
      dxmothet == 7 ~ 19,
      dxmothet == 8 ~ 5,
      dxmothet == 9 ~ 8,
      dxmothet == 10 ~ 12,
      dxmothet == 11 ~ 17,
      dxmothet == 12 ~ 7,
      dxmothet == 13 ~ 1,
      str_detect(dxmothsp, "Multiple Systems Atrophy") ~ 3,
      dxmothet == 14 ~ 30,
      dxmothet == "2|14" ~ 30,
      dxmothet == "5|7|14" ~ 27,
      dxmothet == "7|14" ~ 19,
      dxmothet == "9" ~ 8,
      dxmothet == "9|14" ~ 8,
      dxmothet == "12" ~ 7,
      dxmothet == "7:9" ~ 19,
      dxmothet == "2:14" ~ 18,
      dxmothet == "9:14" ~ 8,
      dxodes == 1 ~ 7, 
      dxodes == 2 ~ 18, 
      dxodes == 3 ~ 11, 
      dxodes == 4 ~ 4, 
      dxodes == 5 ~ 26, 
      dxodes == 6 ~ 14, 
      dxodes == 7 ~ 19, 
      dxodes == 8 ~ 5, 
      dxodes == 9 ~ 8, 
      dxodes == 10 ~ 12, 
      dxodes == 11 ~ 17, 
      dxodes == 12 ~ 7, 
      dxodes == 13 ~ 1, 
      str_detect(dxoothsp, "Lewy") ~ 2,
      dxodes == 14 ~ 30, 
    )
  ) %>%
  select(rid, ptid, phase, viscode, viscode2, examdate, naccudsd, naccetpr)

write.csv(adni_novel, "adni_novel.csv")

## Diagnosis
merge.raw <- read_csv("~/Dropbox/Research/Data/ADNI/ADNIMERGE_25Sep2023.csv") %>%
  janitor::clean_names()


######################################
# Make the tables with this new diagnosis variable
# Case 1. Caide computed with no cholesterol
# Case 2. Caide computed with cholesterol without the chol missing values
# Case 3. Caide computed with imputed cholesterol
######################################

adni_nacc_raw <- 
  read_csv("~/Desktop/adni_nacc/results 11-03-2023/adni_nacc.csv") %>%
  janitor::clean_names()

janitor::get_dupes(adni_nacc_raw) # no dupes

# Make new variables for the tables of Case 1, 2, 3:
table_vars <- adni_nacc_raw %>%
  mutate(
    Gender = gender,
    Diag = case_when(
      diag == 1 ~ "ADRD, MCI",
      diag == 0 ~ "Cog Normal",
      TRUE ~ NA_character_
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
    Hypercholesterolemia3 = case_when(
      (hypchol3 == 1) ~ "Yes",
      (hypchol3 == 0) ~ "No",
      is.na(hypchol3) ~ NA
    ),
    Hypercholesterolemia2 = case_when(
      (hypchol2 == 1) ~ "Yes",
      (hypchol2 == 0) ~ "No",
      is.na(hypchol2) ~ NA
    )
  ) # end of table_vars

###############################
# Make table for case 2:
###############################

adni_nacc_case2 <- table_vars  %>%
  mutate(
    caide_apoe2 = fct_relevel(caide_apoe2, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe2 = fct_relevel(mcaide_apoe2, "Low, e2+", "Low, e3/e3", "Low, e4+",
                               "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                               "High, e2+", "High, e3/e3", "High, e4+")
  ) %>%
  dplyr::select(race, Diag,
                Hypercholesterolemia2, 
                cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                caide2, apoe, caide_apoe2, mcaide2, mcaide_apoe2 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", Diag ~ "categorical",
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
                 race ~ "Race", Diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide2 ~ "CAIDE_w_chol_noNA", mcaide2 ~ "mCAIDE_w_chol_noNA",
                 apoe ~ "APOE", caide_apoe2 ~ "CAIDE x APOE_w_chol_noNA",
                 mcaide_apoe2 ~ "mCAIDE x APOE_w_chol_noNA",
                 Hypercholesterolemia2 ~ "Hypercholesterolemia_w_chol_noNA"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**Caide with no NAs in Cholesterol**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "~/Desktop/adni_nacc/table_adni_nacc_2.docx", 
    vwidth = 648, vheight = 288
  )

###############################
# Make table for case 3:
###############################

adni_nacc_case3 <- table_vars  %>%
  mutate(
    caide_apoe3 = fct_relevel(caide_apoe3, "Low, e2+", "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                              "High, e2+", "High, e3/e3", "High, e4+"),
    mcaide_apoe3 = fct_relevel(mcaide_apoe3, "Low, e2+", "Low, e3/e3", "Low, e4+",
                               "Mid, e2+", "Mid, e3/e3", "Mid, e4+",
                               "High, e2+", "High, e3/e3", "High, e4+")
  ) %>%
  dplyr::select(race, Diag,
                Hypercholesterolemia3, 
                cohort,
                Gender, Age, Education, Obesity, Hypertension, 
                caide3, apoe, caide_apoe3, mcaide3, mcaide_apoe3 
  ) %>%
  tbl_summary(
    by = race,
    type = list(cohort ~ "categorical", Diag ~ "categorical",
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
                 race ~ "Race", Diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide3 ~ "CAIDE_w_imputed_chol", 
                 mcaide3 ~ "mCAIDE_w_imputed_chol",
                 apoe ~ "APOE", caide_apoe3 ~ "CAIDE x APOE_w_imputed_chol",
                 mcaide_apoe3 ~ "mCAIDE x APOE_w_imputed_chol",
                 Hypercholesterolemia3 ~ "Hypercholesterolemia_w_imputed_chol"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**Caide with imputed Cholesterol**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "~/Desktop/adni_nacc/table_adni_nacc_3.docx", 
    vwidth = 648, vheight = 288
  )

###############################
# Make table for case 1:
###############################

adni_nacc_case1 <- table_vars  %>%
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
                #Hypercholesterolemia1 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", Diag ~ "Diagnosis", cohort ~ "Cohort",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                 caide1 ~ "CAIDE_w_no_chol", 
                 mcaide1 ~ "mCAIDE_w_no_chol",
                 apoe ~ "APOE", caide_apoe1 ~ "CAIDE x APOE_w_no_chol",
                 mcaide_apoe1 ~ "mCAIDE x APOE_w_no_chol"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**Caide with no Cholesterol**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "~/Desktop/adni_nacc/table_adni_nacc_1.docx", 
    vwidth = 648, vheight = 288
  )

########################################################
# Regression models with the new diagnosis variable
########################################################

# Set reference
adni_nacc_raw$apoe <- factor(adni_nacc_raw$apoe, ordered = FALSE )
adni_nacc_raw$apoe = relevel(adni_nacc_raw$apoe, ref = "e3/e3")

# Regression Stratified by Race 
adni_nacc_novel_nhw <-  dplyr::filter(adni_nacc_raw, race == "NHW")
adni_nacc_novel_black <-  dplyr::filter(adni_nacc_raw, race == "Black")
adni_nacc_novel_asian <-  dplyr::filter(adni_nacc_raw, race == "Asian")
adni_nacc_novel_his <-  dplyr::filter(adni_nacc_raw, race == "Hispanic")

# 1. logistic model: diagnosis ~ caide + apoe for NHW
model_caide1_apoe_nhw <- glm(diag ~ z_caide1 + apoe, 
                             data=adni_nacc_novel_nhw, family = binomial)
model_caide2_apoe_nhw <- glm(diag ~ z_caide2 + apoe, 
                             data=adni_nacc_novel_nhw, family = binomial)
model_caide3_apoe_nhw <- glm(diag ~ z_caide3 + apoe, 
                             data=adni_nacc_novel_nhw, family = binomial)

# Glance for #1 model for NHW
add_glance_source_note(tbl_regression(model_caide1_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_nhw.png")

# 1. logistic model: diagnosis ~ caide + apoe for Black
model_caide1_apoe_black <- glm(diag ~ z_caide1 + apoe, 
                               data=adni_nacc_novel_black, family = binomial)
model_caide2_apoe_black <- glm(diag ~ z_caide2 + apoe, 
                               data=adni_nacc_novel_black, family = binomial)
model_caide3_apoe_black <- glm(diag ~ z_caide3 + apoe, 
                               data=adni_nacc_novel_black, family = binomial)

# Glance for #1 model for Black
add_glance_source_note(tbl_regression(model_caide1_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_black.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_black.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_black.png")

# 1. logistic model: diagnosis ~ caide + apoe for Asian
model_caide1_apoe_asian <- glm(diag ~ z_caide1 + apoe, 
                               data=adni_nacc_novel_asian, family = binomial)
model_caide2_apoe_asian <- glm(diag ~ z_caide2 + apoe, 
                               data=adni_nacc_novel_asian, family = binomial)
model_caide3_apoe_asian <- glm(diag ~ z_caide3 + apoe, 
                               data=adni_nacc_novel_asian, family = binomial)

# Glance for #1 model for Asian
add_glance_source_note(tbl_regression(model_caide1_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_asian.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_asian.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_asian.png")

# 1. logistic model: diagnosis ~ caide + apoe for Hispanic
model_caide1_apoe_his <- glm(diag ~ z_caide1 + apoe, 
                             data=adni_nacc_novel_his, family = binomial)
model_caide2_apoe_his <- glm(diag ~ z_caide2 + apoe, 
                             data=adni_nacc_novel_his, family = binomial)
model_caide3_apoe_his <- glm(diag ~ z_caide3 + apoe, 
                             data=adni_nacc_novel_his, family = binomial)

# Glance for #1 model for Hispanic
add_glance_source_note(tbl_regression(model_caide1_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_his.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_his.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_his.png")

# 2. logistic model: diagnosis ~ mcaide + apoe for NHW
model_mcaide1_apoe_nhw <- glm(diag ~ z_mcaide1 + apoe, 
                              data=adni_nacc_novel_nhw, family = binomial)
model_mcaide2_apoe_nhw <- glm(diag ~ z_mcaide2 + apoe, 
                              data=adni_nacc_novel_nhw, family = binomial)
model_mcaide3_apoe_nhw <- glm(diag ~ z_mcaide3 + apoe, 
                              data=adni_nacc_novel_nhw, family = binomial)

# Glance for #2 model for NHW
add_glance_source_note(tbl_regression(model_mcaide1_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_nhw.png")

# 2. logistic model: diagnosis ~ mcaide + apoe for Black
model_mcaide1_apoe_black <- glm(diag ~ z_mcaide1 + apoe, 
                                data=adni_nacc_novel_black, 
                                family = binomial)
model_mcaide2_apoe_black <- glm(diag ~ z_mcaide2 + apoe, 
                                data=adni_nacc_novel_black, 
                                family = binomial)
model_mcaide3_apoe_black <- glm(diag ~ z_mcaide3 + apoe, 
                                data=adni_nacc_novel_black, 
                                family = binomial)

# Glance for #2 model for Black
add_glance_source_note(tbl_regression(model_mcaide1_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_black.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_black.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_black.png")

# 2. logistic model: diagnosis ~ mcaide + apoe for Asian
model_mcaide1_apoe_asian <- glm(diag ~ z_mcaide1 + apoe, 
                                data=adni_nacc_novel_asian, 
                                family = binomial)
model_mcaide2_apoe_asian <- glm(diag ~ z_mcaide2 + apoe, 
                                data=adni_nacc_novel_asian, 
                                family = binomial)
model_mcaide3_apoe_asian <- glm(diag ~ z_mcaide3 + apoe, 
                                data=adni_nacc_novel_asian, 
                                family = binomial)

# Glance for #2 model for Asian
add_glance_source_note(tbl_regression(model_mcaide1_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_asian.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_asian.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_asian.png")

# 2. logistic model: diagnosis ~ mcaide + apoe for Hispanic
model_mcaide1_apoe_his <- glm(diag ~ z_mcaide1 + apoe, data=adni_nacc_novel_his, 
                              family = binomial)
model_mcaide2_apoe_his <- glm(diag ~ z_mcaide2 + apoe, data=adni_nacc_novel_his, 
                              family = binomial)
model_mcaide3_apoe_his <- glm(diag ~ z_mcaide3 + apoe, data=adni_nacc_novel_his, 
                              family = binomial)

# Glance for #2 model for Hispanic
add_glance_source_note(tbl_regression(model_mcaide1_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_his.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_his.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_his.png")

# 3 logistic model: diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
#                            apoe + cohort + chol for NHW
model_all3_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol3, 
                      data=adni_nacc_novel_nhw, family = binomial)

model_all2_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol2, 
                      data=adni_nacc_novel_nhw, family = binomial)

model_all1_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort, 
                      data=adni_nacc_novel_nhw, family = binomial)

# glance for all: NHW
add_glance_source_note(tbl_regression(model_all1_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort 
                        (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_nhw.png")

add_glance_source_note(tbl_regression(model_all2_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol2 
                        (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_nhw.png")

add_glance_source_note(tbl_regression(model_all3_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol3 
                        (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_nhw.png")

# 3 logistic model: diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
#                            apoe + cohort + chol for Black
model_all3_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol3, 
                        data=adni_nacc_novel_black, family = binomial)

model_all2_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol2, 
                        data=adni_nacc_novel_black, family = binomial)

model_all1_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                        data=adni_nacc_novel_black, family = binomial)

# glance for all: Black
add_glance_source_note(tbl_regression(model_all1_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_black.png")

add_glance_source_note(tbl_regression(model_all2_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_black.png")

add_glance_source_note(tbl_regression(model_all3_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_black.png")

# 3 logistic model: diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
#                            apoe + cohort + chol for Asian

model_all3_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol3, 
                        data=adni_nacc_novel_asian, family = binomial)

model_all2_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol2, 
                        data=adni_nacc_novel_asian, family = binomial)

model_all1_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                        data=adni_nacc_novel_asian, family = binomial)

# glance for all: Asian
add_glance_source_note(tbl_regression(model_all1_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_asian.png")

add_glance_source_note(tbl_regression(model_all2_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_asian.png")

add_glance_source_note(tbl_regression(model_all3_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_asian.png")

# 3 logistic model: diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
#                            apoe + cohort + chol for Hispanic
model_all3_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                        gender + apoe + cohort + hypchol3, 
                      data=adni_nacc_novel_his, family = binomial)

model_all2_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                        gender + apoe + cohort + hypchol2, 
                      data=adni_nacc_novel_his, family = binomial)

model_all1_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                        gender + apoe + cohort, 
                      data=adni_nacc_novel_his, family = binomial)

# glance for all: Hispanic
add_glance_source_note(tbl_regression(model_all1_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_his.png")


add_glance_source_note(tbl_regression(model_all2_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_his.png")

add_glance_source_note(tbl_regression(model_all3_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_his.png")




