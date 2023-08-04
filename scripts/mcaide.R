# Filename: mcaide_last_v2.R
# calculation of caide & mcaide

setwd("~/Desktop/Chi_data3/data")
library("tidyverse")
library("gtsummary")
library("missForest")
# library("dplyr")
library("glue")
library("magrittr")
library("MASS")

chidata_orig <- read_csv("chi_data1.csv", guess_max = 10000) %>% 
  janitor::clean_names()
# procount(chidata_orig, diag_x)

chi_orig <- subset(chidata_orig, select = -c(apoe, cohort))

# complete diagnosis:
chi_diag <- chi_orig %>%
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
  ) 

chi_diag %>% count(diagnosisvisit_1)
# CI: 273, N: 4080

# complete apoe
apoe_wide <- chi_diag %>%
  mutate(
    apoesummary = case_when(
      !is.na(apoevisit_2) ~ apoevisit_2,
      is.na(apoevisit_2) & apoevisit_1 == "E4 carrier" ~ "44",
      TRUE ~ NA_character_ 
    )
  ) 

# fix typing error
apoe_wide[apoe_wide == "4Y"] <- "44"

# attach apoe data to lv data

# pivot data from wide to long (& then slice)
chi_long <- apoe_wide %>%
  dplyr::select(subid, apoesummary, starts_with('cohort'),starts_with('time'), starts_with('cog'), 
                starts_with('age'), starts_with('diag')) %>%
  dplyr::select(!dplyr::contains('last')) %>% 
  pivot_longer(
    cols = !c(subid, apoesummary),
    names_to = c(".value", "visit"),
    names_pattern = "(.*)visit_(.*)"
  ) %>%
  dplyr::select(-ageatdeath) 

chi_long %>% count(apoesummary)

# filter to last visit patient attended (slice):
chi_lv <- chi_long %>% 
  group_by(subid) %>%
  slice_max(time) %>%
  ungroup() %>%
  magrittr::set_colnames(paste0(colnames(.), "_lv"))

# check diagnosis:
chi_lv_new_diag <-  dplyr::filter(chi_lv, cohort_lv == 'New')
chi_lv_new_diag %>% count(diagnosis_lv)
# diagnosis_lv: CI: 240, D: 64, N: 1591

# append last visit to wide 
# check apoe:
chi_lv_new <-  dplyr::filter(chi_lv, cohort_lv == 'New')
chi_lv_new %>% count(apoesummary_lv) # apoesum: 887, NA: 1008

chi_diag_new <- dplyr::filter(chi_diag, cohortvisit_1 == 'New')
chi_diag_new %>% count(apoevisit_2) #apoe: 833, 4Y: 1, NA: 1061

# left join chi_lv (last visit) and chi_diag (wide):
names(chi_lv)[names(chi_lv) == 'subid_lv'] <- 'subid'
lv_to_wide <- merge(x = chi_lv, y = chi_diag, by = 'subid', 
                    all.x=TRUE)

# restrict to 'New' cohort:
chi_combined_new <-  dplyr::filter(lv_to_wide, cohort_lv == 'New')
chi_combined_new %>% count(apoesummary_lv) # apoesum: 887, NA: 1008
chi_combined_new %>% count(diagnosis_lv) # CI: 240, D: 64, N: 1591

tail(chi_combined_new) # looks good
# dplyr::select(chi_long_max_time, subid, visit, time, diagnosis) %>% 
# print(n=100)

# drop variables that are not baseline or lv, 
#    retain any variable for caide at baseline
####
# select columns:
selected_chi <- dplyr::select(chi_combined_new, c('subid', 'apoesummary_lv', 
                                                  'visit_lv', 'cohort_lv',
                                                  'time_lv', 'cogfull_lv',
                                                  'cog_cat_lv', 'age_lv',
                                                  'diag_lv', 'diagnosis_lv',
                                                  'bmivisit_1', 'systbpvisit_1',
                                                  'sexvisit_1', 'gradevisit_1',
                                                  'agevisit_1', 
                                                  'physical_act'))

str(selected_chi) 
selected_chi %>% count(apoesummary_lv) # apoesummary: 887
tail(selected_chi)

# prepare data for imputing missing values:
chi_NoNAs4apoe <- subset(selected_chi, apoesummary_lv != "")
chi_NoNAs4apoe %>% count(apoesummary_lv) # looks good
tail(chi_NoNAs4apoe)

# convert categorical to factors:
chi_factors <- chi_NoNAs4apoe %>%
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

# prepare missing values calculation
chi_selected <- subset(chi_factors, select = c( 
                                                    sex_factor,
                                                    gradevisit_1,
                                                    agevisit_1, 
                                                    bmivisit_1,
                                                    systbpvisit_1,
                                                    phys_act_factor
                                                  ))
tail(chi_selected) 
chi_selected %>% count(apoe_factor)
chi_selected %>% count(diagnosis_factor) # diagnosis_factor not present

chi_selected_df <- as.data.frame(chi_selected) 
chi_factor_only <- subset(chi_factors, select = -c(
                                                       sex_factor,
                                                       gradevisit_1,
                                                       agevisit_1,
                                                       bmivisit_1,
                                                       systbpvisit_1,
                                                       phys_act_factor
                                                      ))
 
# imputation, using only variables for caide
no_missing_values <- missForest(chi_selected_df)
chi_nmv <- no_missing_values$ximp
tail(chi_nmv)

# join chi_nmv and chi_factor_only:
chi_nmv_tib <- tibble(chi_nmv)
tibble(chi_nmv_tib)
chi_factor_tib <- tibble(chi_factor_only)
tibble(chi_factor_tib)

chi_nmv_cb <- cbind(chi_nmv_tib, chi_factor_tib)

# clean up data (fix i_names):
chidata <- chi_nmv_cb %>%
  rename("gender" = "sexvisit_1") %>%
  rename("i_age" = "agevisit_1") %>%
  rename("i_education" = "gradevisit_1") %>%
  rename("i_bmi" = "bmivisit_1") %>%
  rename("i_systbp" = "systbpvisit_1") %>%
  rename("apoe" = "apoe_factor") %>%
  rename("diagnosis" = "diagnosis_lv") %>%
  rename("cognitive" = "cogfull_lv")

chidata %>% count(apoe)
chidata %>% count(diagnosis)
# calculation of caide and mcaide:
# 1 compute caide_sex, caide_age, etc
# 2 compute caide, caide_cat, etc 
# 3 convert apoe alleles from "22", "23", etc to "e2+", "e3/e3", etc
# 4 compute caide x apoe, etc
# 5 compute the above for mcaide

# compute caide_sex, caide_age, etc:
caide_1 <- chidata %>%
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
  )

# compute caide, caide_cat, etc; convert apoe alleles from "22", "23", etc to "e2+", "e3/e3", etc
caide_2 <- caide_1 %>%
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
            between(caide, 0, (mean(caide, na.rm = T) - sd(caide, na.rm = T))) ~ 'Low',
            between(caide, (mean(caide, na.rm = T) - sd(caide, na.rm = T)), (mean(caide, na.rm = T) + sd(caide, na.rm = T))) ~ 'Mid',
            between(caide, (mean(caide, na.rm = T) + sd(caide, na.rm = T)), 14) ~ 'High',
            TRUE ~ NA_character_
          )
  ) %>%
  mutate(
    apoe_ = case_when(
      apoe == 22 ~ 'e2+',
      apoe == 23 ~ 'e2+',
      apoe == 24 ~ 'e4+',
      apoe == 33 ~ 'e3/e3',
      apoe == 34 ~ 'e4+',
      apoe == 44 ~ 'e4+',
      TRUE ~ 'NA'
    )
  ) %>%
  mutate(
    caide_cat = fct_relevel(caide_cat, 'Low', "Mid",  'High'),
    caide_apoe = glue("{caide_cat}, {apoe_}"),
    caide_apoe = fct_relevel(caide_apoe, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                             "Low, e4+",
                             "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                             "High, e4+")
  ) %>%
  mutate_at(
    vars(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
         caide_phy_inac), as_factor
  ) %>%
  dplyr::select(subid, starts_with('caide'), z_caide, caide_cat, apoe_, 
                caide_apoe, i_age, i_education, gender, i_bmi, i_systbp, 
                physical_act, caide_phy_inac, diagnosis, cognitive)

write_csv(x = caide_2, "caide_v2_last.csv")

# calculation of mcaide:

mcaide1 <- chidata %>%
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
      #physical_act == 0 ~ 1,
      #physical_act == 1 ~ 0,
      physical_act == 0 ~ 1,
      physical_act == 1 ~ 0,
      TRUE ~ NA_real_
    ),
  )

mcaide_2 <- mcaide1 %>%
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
    apoe_ = case_when(
      apoe == 22 ~ 'e2+',
      apoe == 23 ~ 'e2+',
      apoe == 24 ~ 'e4+',
      apoe == 33 ~ 'e3/e3',
      apoe == 34 ~ 'e4+',
      apoe == 44 ~ 'e4+',
      TRUE ~ 'NA'
    ),
    mcaide_cat = fct_relevel(mcaide_cat, 'Low', "Mid",  'High'),
    mcaide_apoe = glue("{mcaide_cat}, {apoe_}"),
    mcaide_apoe = fct_relevel(mcaide_apoe, "Mid, e3/e3", "Low, e2+", 
                              "Low, e3/e3", "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", 
                              "High, e3/e3", "High, e4+")
  ) %>%
  mutate_at(
    vars(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp, 
         mcaide_phy_inac), as_factor
  ) %>%
  dplyr::select(subid, starts_with('mcaide'), z_mcaide, mcaide_cat, apoe_, 
                mcaide_apoe, i_age, i_education, gender, i_bmi, i_systbp, physical_act, 
                mcaide_phy_inac, diagnosis, cognitive)

write_csv(x = mcaide_2, "mcaide_v2_last.csv")

# make table --> see diag.txt

#regression models:
library("MASS")

data <- read_csv("mcaide_v2_last.csv")
data_caide <- read_csv("caide_v2_last.csv")
data$diagnosis <- as.factor(data$diagnosis)

chi_data_mcaide <- data %>%
  mutate(
    diag_ad = case_when(
      diagnosis == 'N' ~ 0,
      diagnosis == 'CI' ~ 1,
      diagnosis == 'D' ~ 1,
    )
  )

chi_data_mcaide$apoe_ <- as.factor(chi_data_mcaide$apoe_)
chi_data_mcaide$apoe_ = relevel(chi_data_mcaide$apoe_, ref = "e3/e3")
model_mcaide_log <- glm(data=chi_data_mcaide, diag_ad ~ apoe_ + mcaide,
                      family="binomial")
summary(model_mcaide_log)

data_caide$diagnosis <- as.factor(data_caide$diagnosis)

chi_data_caide <- data_caide %>%
  mutate(
    diag = case_when(
      diagnosis == 'N' ~ 0,
      diagnosis == 'CI' ~ 1,
      diagnosis == 'D' ~ 1,
    )
  )

chi_data_caide$apoe_ <- as.factor(chi_data_caide$apoe_)
chi_data_caide$apoe_ = relevel(chi_data_caide$apoe_, ref = "e3/e3")
model_caide_log <- glm(data=chi_data_caide, diag ~ apoe_ + caide, 
                    family="binomial")
summary(modelcaide_log)

# linear regression with cogfull the outcome, mcaide
chi_data_mcaide$apoe_ <- factor(chi_data_mcaide$apoe_ , 
                                   ordered = FALSE )
chi_data_mcaide$apoe_ = relevel(chi_data_mcaide$apoe_, 
                                   ref = "e3/e3")
model_mcaide_cogfull <- lm(data=chi_data_mcaide, cognitive ~ apoe_ + 
                            mcaide)
summary(model_mcaide_cogfull)

library("broom")
add_glance_source_note(tbl_regression(model_caide_log)) # gives error
broom::glance(model_mcaide_log)

# linear regression with cogfull the outcome, caide
chi_data_caide$apoe_ <- factor(chi_data_caide$apoe_ , 
                         ordered = FALSE )
chi_data_caide$apoe_ = relevel(chi_data_caide$apoe_, 
                         ref = "e3/e3")
model_caide_cogfull <- lm(data=chi_data_caide, cognitive ~ apoe_ + 
                             caide)
summary(model_caide_cogfull)
