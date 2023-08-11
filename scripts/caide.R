# Filename: caide.R
# Date: 08-11-2023
# calculation of caide & mcaide after munging data

library("tidyverse")
library("gtsummary")
library("missForest")
library("dplyr")
library("glue")
library("magrittr")
library("MASS")
library("readxl")

data <- read_excel("CAIDE_Risk_Score_and_APOE_project_with_Shea.xlsx", 
                   sheet = "Sheet1", guess_max = 10000) %>% 
                   janitor::clean_names()

# complete diagnosis:
iidp_wide <- subset(data, select = -c(apoe, cohort)) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "\\_x", "visit_1")) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "\\_y", "visit_2")) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "visit_1visit_1", "visit_3")) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "visit_2visit_2", "visit_4")) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "visit_3visit_1", "visit_5")) %>%
  magrittr::set_colnames(str_replace_all(colnames(.), "visit_4visit_2", "visit_6")) %>%
  mutate(
    diagnosisvisit_1 = case_when(
      !is.na(timevisit_1) & is.na(diagvisit_1) ~ "N",
      !is.na(timevisit_1) & !is.na(diagvisit_1) ~ diagvisit_1,
      TRUE ~ NA_character_
    ),
    diagnosisvisit_2 = case_when(
      !is.na(timevisit_2) & is.na(diagvisit_2) ~ "N",
      !is.na(timevisit_2) & !is.na(diagvisit_2) ~ diagvisit_2,
      TRUE ~ NA_character_
    ),
    diagnosisvisit_3 = case_when(
      !is.na(timevisit_3) & is.na(diagvisit_3) ~ "N",
      !is.na(timevisit_3) & !is.na(diagvisit_3) ~ diagvisit_3,
      TRUE ~ NA_character_
    ),
    diagnosisvisit_4 = case_when(
      !is.na(timevisit_4) & is.na(diagvisit_4) ~ "N",
      !is.na(timevisit_4) & !is.na(diagvisit_4) ~ diagvisit_4,
      TRUE ~ NA_character_
    ),
    diagnosisvisit_5 = case_when(
      !is.na(timevisit_5) & is.na(diagvisit_5) ~ "N",
      !is.na(timevisit_5) & !is.na(diagvisit_5) ~ diagvisit_5,
      TRUE ~ NA_character_
    ),
    diagnosisvisit_6 = case_when(
      !is.na(timevisit_6) & is.na(diagvisit_6) ~ "N",
      !is.na(timevisit_6) & !is.na(diagvisit_6) ~ diagvisit_6,
      TRUE ~ NA_character_
    ) 
  ) %>%
# complete apoe
  mutate(
    apoesummary = case_when(
      !is.na(apoevisit_2) ~ apoevisit_2,
      is.na(apoevisit_2) & apoevisit_1 == "E4 carrier" ~ "44",
      TRUE ~ NA_character_ 
    )
  ) 

# fix data error
iidp_wide[iidp_wide == "4Y"] <- "44"

# attach apoe data to lv data
# pivot data from wide to long
iidp_long <- iidp_wide %>%
  dplyr::select(subid, apoesummary, starts_with('cohort'),starts_with('time'), starts_with('cog'), 
                starts_with('age'), starts_with('diag')) %>%
  dplyr::select(!dplyr::contains('last')) %>% 
  pivot_longer(
    cols = !c(subid, apoesummary),
    names_to = c(".value", "visit"),
    names_pattern = "(.*)visit_(.*)"
  ) %>%
  dplyr::select(-ageatdeath) 

iidp_long %>% count(apoesummary)

# filter to last visit patient attended (slice):
iidp_lv <- iidp_long %>% 
  group_by(subid) %>%
  slice_max(time) %>%
  ungroup() %>%
  magrittr::set_colnames(paste0(colnames(.), "_lv"))

# check diagnosis in 'New' cohort in last visit file:
iidp_lv_New_diag <-  dplyr::filter(iidp_lv, cohort_lv == 'New')
iidp_lv_New_diag %>% count(diagnosis_lv)
# diagnosis_lv: CI: 240, D: 64, N: 1591

# check apoe in 'New' cohort in last visit file:
iidp_lv_New_apoe <-  dplyr::filter(iidp_lv, cohort_lv == 'New')
iidp_lv_New_apoe %>% count(apoesummary_lv) # apoesum: 887, NA: 1008

# append last visit to wide:
# left join iidp_lv (last visit) and iidp_wide (wide):
names(iidp_lv)[names(iidp_lv) == 'subid_lv'] <- 'subid'

lv_to_wide <- merge(x = iidp_lv, y = iidp_wide, by = 'subid', 
                    all.x=TRUE)

# restrict to 'New' cohort:
iidp_merged_New <-  dplyr::filter(lv_to_wide, cohort_lv == 'New')

# check apoe and diagnosis in left joined file:
iidp_merged_New %>% count(apoesummary_lv) # apoesum: 887, NA: 1008
iidp_merged_New %>% count(diagnosis_lv) # CI: 240, D: 64, N: 1591
tail(iidp_merged_New) # looks good

# drop variables that are not baseline or lv, 
#    retain any variable for caide at baseline

# select columns:
selected_iidp <- dplyr::select(iidp_merged_New, c('subid', 'apoesummary_lv', 
                                                  'visit_lv', 'cohort_lv',
                                                  'time_lv', 'cogfull_lv',
                                                  'cog_cat_lv', 'age_lv',
                                                  'diag_lv', 'diagnosis_lv',
                                                  'bmivisit_1', 'systbpvisit_1',
                                                  'sexvisit_1', 'gradevisit_1',
                                                  'agevisit_1', 
                                                  'physical_act'))

# check diagnosis in selected file:
selected_iidp %>% count(diagnosis_lv) # diagnosis: CI: 240, D: 64, N:1591

# prepare data for imputing missing values:
iidp_NoNAs4apoe <- subset(selected_iidp, apoesummary_lv != "")

# check diagnosis in 'NoNAs' file: Lost 50% of the data!
iidp_NoNAs4apoe %>% count(diagnosis_lv) # diagnosis: CI: 109, D: 37, N: 741

# convert categorical to factors:
iidp_factors <- iidp_NoNAs4apoe %>%
  mutate(
         apoe_factor = as.factor(apoesummary_lv),
	       visit_factor = as.factor(visit_lv),
	       cohort_factor = as.factor(cohort_lv),
	       cog_cat_factor = as.factor(cog_cat_lv),
         sex_factor = as.factor(sexvisit_1),
	       diag_factor = as.factor(diag_lv),
         diagnosis_factor = as.factor(diagnosis_lv),
         phys_act_factor = as.factor(physical_act)
  )

# split up file into (1) imputation file and (2) other values file:
iidp_imputation <- subset(iidp_factors, select = c( 
                                                    sex_factor,
                                                    gradevisit_1,
                                                    agevisit_1, 
                                                    bmivisit_1,
                                                    systbpvisit_1,
                                                    phys_act_factor
                                                  ))

iidp_factor_only <- subset(iidp_factors, select = -c(
                                                       sex_factor,
                                                       gradevisit_1,
                                                       agevisit_1,
                                                       bmivisit_1,
                                                       systbpvisit_1,
                                                       phys_act_factor
                                                      ))
 
# imputation, using only variables for caide
iidp_imputation_df <- as.data.frame(iidp_imputation)
no_missing_values <- missForest(iidp_imputation_df)
iidp_nmv <- no_missing_values$ximp

# check number of non-null rows in tibbles:
iidp_nmv_tib <- tibble(iidp_nmv)
iidp_factor_tib <- tibble(iidp_factor_only)
str(iidp_nmv_tib)
str(iidp_factor_tib)

# join chi_nmv and chi_factor_only by column bind:
iidp_nmv_cb <- cbind(iidp_nmv_tib, iidp_factor_tib)

# fix i_names:
iidpdata <- iidp_nmv_cb %>%
  rename("gender" = "sexvisit_1") %>%
  rename("i_age" = "agevisit_1") %>%
  rename("i_education" = "gradevisit_1") %>%
  rename("i_bmi" = "bmivisit_1") %>%
  rename("i_systbp" = "systbpvisit_1") %>%
  rename("apoe_f" = "apoe_factor") %>%
  rename("diagnosis" = "diagnosis_lv") %>%
  rename("cognitive" = "cogfull_lv")

iidpdata %>% count(apoe_f) # ok
iidpdata %>% count(diagnosis) # CI: 109, D: 37, N: 741

# compute caide_sex, caide_age, etc (two parts):
caide_calc <- iidpdata %>%
  mutate(
    caide_age = case_when(
      i_age <= 46 ~ 0,
      i_age >= 47 & i_age < 53 ~ 3,
      i_age >= 53 ~ 4,
      TRUE ~ NA_real_
    ),
    caide_educ = case_when(
      i_education <= 6 ~ 3,
      i_education >= 7 & i_education <= 9 ~ 2,
      i_education > 9 ~ 0,
      TRUE ~ NA_real_
    ),
    caide_sex = case_when(
      gender == "M" ~ 1,
      gender == "F" ~ 0,
      TRUE ~ NA_real_
    ),
    caide_obesity = case_when(
      i_bmi < 30 ~ 0,
      i_bmi >= 30 ~ 2,
      TRUE ~ NA_real_
    ),
    caide_sbp = case_when(
      i_systbp < 140 ~ 0,
      i_systbp >= 140 ~ 2,
      TRUE ~ NA_real_
    ),
    caide_phy_inac = case_when(
      physical_act == 0 ~ 1,
      physical_act == 1 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    caide = sum(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
                caide_phy_inac, na.rm = F),
    caide_missing = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                        is.na(i_bmi), is.na(i_systbp), is.na(physical_act))
  ) %>%
  ungroup() %>%
  mutate(
    z_caide = scale(caide)[,1],
    caide_cat = case_when(
            between(caide, 0, (mean(caide, na.rm = T) - sd(caide, na.rm = T))) 
            ~ 'Low',
            between(caide, (mean(caide, na.rm = T) - sd(caide, na.rm = T)), 
                    (mean(caide, na.rm = T) + sd(caide, na.rm = T))) ~ 'Mid',
            between(caide, (mean(caide, na.rm = T) + sd(caide, na.rm = T)), 14) ~ 'High',
            TRUE ~ NA_character_
          )
  ) %>%
  mutate(
    apoe = case_when(
      apoe_f == 22 ~ 'e2+',
      apoe_f == 23 ~ 'e2+',
      apoe_f == 24 ~ 'e4+',
      apoe_f == 33 ~ 'e3/e3',
      apoe_f == 34 ~ 'e4+',
      apoe_f == 44 ~ 'e4+',
      TRUE ~ 'NA'
    )
  ) %>%
  mutate(
    caide_cat = fct_relevel(caide_cat, 'Low', "Mid",  'High'),
    caide_apoe = glue("{caide_cat}, {apoe}"),
    caide_apoe = fct_relevel(caide_apoe, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                             "Low, e4+",
                             "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                             "High, e4+")
  ) %>%
  mutate_at(
    vars(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
         caide_phy_inac), as_factor
  ) %>%
  dplyr::select(subid, starts_with('caide'), z_caide, caide_cat, apoe, 
                caide_apoe, i_age, i_education, gender, i_bmi, i_systbp, 
                physical_act, caide_phy_inac, diagnosis, cognitive)

write_csv(x = caide_calc, "caide_v2_last.csv")

# calculation of mcaide:

mcaide_calc <- iidpdata %>%
  mutate(
    mcaide_age = case_when(
      i_age <= 64 ~ 0,
      i_age >= 65 & i_age < 73 ~ 1,
      i_age >= 73 ~ 2,
      TRUE ~ NA_real_
    ),
    mcaide_educ = case_when(
      i_education <= 11 ~ 2,
      i_education >= 12 & i_education <= 16 ~ 1,
      i_education > 16 ~ 0,
      TRUE ~ NA_real_
    ),
    mcaide_sex = case_when(
      gender == "M" ~ 1,
      gender == "F" ~ 0,
      TRUE ~ NA_real_
    ),
    mcaide_obesity = case_when(
      i_bmi <= 30 ~ 0,
      i_bmi > 30 ~ 2,
      TRUE ~ NA_real_
    ),
    mcaide_sbp = case_when(
      i_systbp < 140 ~ 0,
      i_systbp >= 140 ~ 2,
      TRUE ~ NA_real_
    ),
    mcaide_phy_inac = case_when(
      physical_act == 0 ~ 1,
      physical_act == 1 ~ 0,
      TRUE ~ NA_real_
    ),
  ) %>%
  rowwise() %>%
  mutate(
    mcaide = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp, mcaide_phy_inac, na.rm = F),
    mcaide_missing = sum(is.na(i_age), is.na(i_education), is.na(gender),
                        is.na(i_bmi), is.na(physical_act), is.na(i_systbp)),
  ) %>%
  ungroup() %>%
  mutate(
    z_mcaide = scale(mcaide)[,1],
    mcaide_cat = case_when(
      between(mcaide, 0, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T))) ~ 'Low',
      between(mcaide, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T)), (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T))) ~ 'Mid',
      between(mcaide, (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T)), 14) ~ 'High',
      TRUE ~ NA_character_
    ) 
  ) %>%
  mutate(
    apoe = case_when(
      apoe_f == 22 ~ 'e2+',
      apoe_f == 23 ~ 'e2+',
      apoe_f == 24 ~ 'e4+',
      apoe_f == 33 ~ 'e3/e3',
      apoe_f == 34 ~ 'e4+',
      apoe_f == 44 ~ 'e4+',
      TRUE ~ 'NA'
    ),
    mcaide_cat = fct_relevel(mcaide_cat, 'Low', "Mid",  'High'),
    mcaide_apoe = glue("{mcaide_cat}, {apoe}"),
    mcaide_apoe = fct_relevel(mcaide_apoe, "Mid, e3/e3", "Low, e2+", 
                              "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", 
                              "High, e3/e3", "High, e4+")
  ) %>%
  mutate_at(
    vars(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp, 
         mcaide_phy_inac), as_factor
  ) %>%
  dplyr::select(subid, starts_with('mcaide'), z_mcaide, mcaide_cat, apoe, 
                mcaide_apoe, i_age, i_education, gender, i_bmi, i_systbp, physical_act, 
                mcaide_phy_inac, diagnosis, cognitive)

write_csv(x = mcaide_calc, "mcaide_v2_last.csv")

# make table --> see table1.R
# do regression --> see regression.R

