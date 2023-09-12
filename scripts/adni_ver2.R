# Filename: adni_ver2.R

# setwd("/wynton/group/andrews/users/evelasquez/CAIDE_APOE/scripts/ADNI")
setwd("/Users/elinorvelasquez/Desktop/ADNI/")
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
library("cowplot")
library("patchwork")

#####################################################
# compute data_final (contains not-to-be-imputed variables)
#####################################################

data <- read_csv("ADNIMERGE_09Aug2023.csv") %>% 
  janitor::clean_names()

data_sc_bl <- dplyr::filter(data, viscode == "sc" | viscode == "bl")
# (missing "race" var:)

#### Add race back in:
data_sc_bl_race <- data_sc_bl %>%
  mutate(
    race = case_when(
      # ptraccat != "Asian" & ptethcat != "Hisp/Latino" & ptraccat != "Black" 
      #  & ptraccat != "White" ~ "Other",
      ptraccat == "More than one" ~ "Other",
      ptraccat == "Unknown" ~ "Other",
      ptraccat == "Am Indian/Alaskan" ~ "Other",
      ptraccat == "Hawaiian/Other PI" ~ "Other",
      ptraccat == "White" & ptethcat != "Hisp/Latino" ~ "Non-Hispanic White",
      ptraccat == "White" & ptethcat == "Hisp/Latino" ~ "Hispanic",
      ptraccat == "Black" ~ "Black",
      ptraccat == "Asian" ~ "Asian",
      TRUE ~ NA_character_
    )
  )
####

data_no_impute_vars <- dplyr::select(data_sc_bl_race, rid, race, dx_bl, 
                                     adas11_bl, m_pac_cdigit_bl)
# no mv for dx_bl
no_mv_for_dx_bl <-  dplyr::filter(data_no_impute_vars, dx_bl != "NA") 
data_final1 <- no_mv_for_dx_bl %>% dplyr::distinct()

# check:
data_final1 %>% count(dx_bl) # yes, all NAs are gone
# need to get rid of "Other"
data_final <- dplyr::filter(data_final1, race != "Other")
data_final %>% count(race) # yes, "Other" is gone?

# For now, ignore the missing values of adas & mpacc
data_final %>% count(!is.na(adas11_bl)) # TRUE 2417, FALSE 13
data_final %>% count(!is.na(m_pac_cdigit_bl)) # TRUE 2427, FALSE 3

###################################################
# Before imputing the caide variables, compute caide_basic & bmi-systbp
###################################################
# 1. Compute caide_basic:

adnimerge_subset <- dplyr::select(data_sc_bl, rid, viscode, age, ptgender, 
                                  pteducat)
caide_basic_fa <- adnimerge_subset %>% 
  mutate(
    viscode_factor = as.factor(viscode),
    sex_factor = as.factor(ptgender)
  )
caide_basic_imp <- dplyr::select(caide_basic_fa, -c(viscode, ptgender))
caide_basic = dplyr::rename(caide_basic_imp, viscode = viscode_factor)

# check for missing values
summary(caide_basic) # age has 4 mv; all others have no mv

# 2. Compute bmi-systbp:

#####################################################
# COMPUTE BMI & SYSBP: Here is "data_bmi_sc2"
#####################################################

# Construct a file with "sc" height, weight, and respective units & systbp. 
# Compute bmi using these values.
vitals <- read_csv("VITALS_14Aug2023.csv") %>% 
  janitor::clean_names()
data_bmi <- dplyr::select(vitals, rid, viscode, viscode2, 
                          vsweight, vswtunit, 
                          vsheight, vshtunit, vsbpsys)

data_bmi_sc <- dplyr::filter(data_bmi, viscode == "sc")
data_bmi_sc %>% count(!is.na(vsbpsys)) # TRUE 2071, FALSE 21
data_bmi_sc2 <- dplyr::filter(data_bmi, viscode2 == "sc") 
data_bmi_sc2 %>% count(!is.na(vsbpsys)) # TRUE 3179, FALSE 32
data_bmi_sc2 %>% count(!is.na(vsweight)) # TRUE 3196, FALSE 15
data_bmi_test <- dplyr::filter(data_bmi, viscode != viscode2)
# we definitely want to use sc and viscode2
# Need to change "-1, -4" entries to "NA":
data_bmi_sc2[data_bmi_sc2 == -1] <- NA_real_
data_bmi_sc2[data_bmi_sc2 == -4] <- NA_real_ 

# Compute BMI
bmi_computed <- data_bmi_sc2 %>%
  mutate(
    weight_in_kg = case_when(vswtunit == 1 ~ vsweight*0.453592, #1 is lbs
                             vswtunit == 2 ~ vsweight,
                             TRUE ~ NA_real_
    ),
    height_in_cm = case_when(vshtunit == 1 ~ vsheight*2.54, #1 is inches
                             vshtunit == 2 ~ vsheight,
                             TRUE ~ NA_real_
    ),
    bmi = weight_in_kg/((height_in_cm)**2)*10000 # From CDC
  )

# 3. compute bmi_systbp:
bmi_systbp_imp <- dplyr::select(bmi_computed, rid, viscode2, bmi, vsbpsys)
bmi_systbp_fa <- bmi_systbp_imp %>% 
  mutate(
    viscode_factor = as.factor(viscode2)
  )
bmi_systbp_fac <- dplyr::select(bmi_systbp_fa, -c(viscode2))
bmi_systbp <- dplyr::rename(bmi_systbp_fac, viscode = viscode_factor)
##### BMI & systbp computed #######################

###################################################
# Impute without cholesterol: Need caide_basic, bmi_systbp
###################################################

caide_basic_join <- dplyr::select(caide_basic, rid, age, pteducat, 
                                  sex_factor)
bmi_systbp_join <- dplyr::select(bmi_systbp, rid, bmi, vsbpsys)

impute_wo_chol_list <- list(caide_basic_join, bmi_systbp_join)
caide_wo_chol_joined <- impute_wo_chol_list %>% reduce(full_join, by='rid')

# check:
summary(caide_basic_join) # only age has missing values --4
summary(bmi_systbp_join) # bmi na's: 44; bpsys na's: 32

left_joined <- caide_basic_join %>% left_join(bmi_systbp_join, 
                             by=c('rid'))
summary(left_joined)
str(left_joined) # 2430 rows # switch to 'left_joined'

left_joined[left_joined == -1] <- NA_real_
left_joined[left_joined == -4] <- NA_real_
summary(left_joined)

# check
str(caide_wo_chol_joined) # rows 3211
summary(caide_wo_chol_joined)

# imputing...
left_joined_df <- as.data.frame(left_joined)
no_missing_values_lj <- missForest(left_joined_df)
lj_nmv <- no_missing_values_lj$ximp
imputed_wo_chol_lj <- tibble(lj_nmv)

# testing ...
summary(imputed_wo_chol_lj) # no mv's
imputed_wo_chol_lj %>% count(!is.na(vsbpsys)) # TRUE 3211 FALSE 0
str(imputed_wo_chol_lj) # rows 2430
                        # contains rid, age, pteducat, sex_factor, bmi, vsbpsys
write_csv(imputed_wo_chol, "imputed_wo_chol.csv")

###########################################
# Caide vars imputed; Will need to add back in chol
###########################################
# Chol dataset constructed:
###################################################################
# COMPUTE CHOLESTEROL ("chol_final" dataset)
###################################################################

chol_data <- read_csv("biospecimen_ADNINIGHTINGALE2_15Aug2023.csv") %>% 
  janitor::clean_names()
# Test if viscode1 == viscode2 in the dataset NIGHTINGALE
chol_data %>% count(viscode == viscode2) # TRUE 928, FALSE 769
chol_test12 <- dplyr::filter(chol_data, viscode == viscode2) 
chol_test_bl <- dplyr::filter(chol_data, viscode2 == "bl") 
chol_test_sc <- dplyr::filter(chol_data, viscode2 == "sc") 
chol_test_sc %>% count(!is.na(total_c)) 
# viscode1: 924 non-zero values for bl
# viscode2: 1693 non-zero values for bl
# viscode1: 1 non-zero value for sc
# viscode2: 1 non-zero value for sc

# test to examine the viscodes:
test_chol <- dplyr::filter(chol_data, viscode != viscode2)
# observed ALL viscode = "v03", while viscode2 = "bl" if viscode1 != viscode2
# Conclusion: use viscode2 for chol
chol_complete <- dplyr::filter(chol_data, viscode2 == "bl" | viscode2 == "sc")

chol_complete %>% count(!is.na(total_c)) # TRUE 1694 FALSE 1
chol_select <- dplyr::select(chol_complete, rid, total_c)
chol_final <- chol_select %>% dplyr::distinct() # remove duplicates
chol_final %>% count(!is.na(total_c)) # TRUE 1694, FALSE 1

####################################################
# Construct APOE now:
####################################################

# APOE doesn't have many values. See if we can fix that. 
apoe_meta_new <- read_csv("ADSP_PHC_ADNI_T1_1.0_MetaData_15Aug2023.csv") %>% 
  janitor::clean_names()
apoe_meta_name = dplyr::rename(apoe_meta_new, viscode = loni_viscode2)
apoe_adsp_selected <- dplyr::select(apoe_meta_name, rid, viscode, apoe4count, 
                                    apoe3count, apoe2count)
apoe_adsp_bl <- dplyr::filter(apoe_adsp_selected, viscode == "bl")
apoe_adsp_bl %>% count(!is.na(apoe4count)) # TRUE 1521, FALSE 836
apoe_adsp_sc <- dplyr::filter(apoe_adsp_selected, viscode == "sc")
apoe_adsp_sc %>% count(!is.na(apoe4count)) # TRUE 0, FALSE 92
# Conclusion: use "bl" only. See if we can merge the two apoe datasets

# New dataset for apoe: (has viscode & rid)
apoe_alt <- read_csv("APOERES_16Aug2023.csv") %>% 
  janitor::clean_names()

apoe_alt %>% count(apgen1, phase)
apoe_alt %>% count(!is.na(apgen2)) # TRUE 2556, FALSE 0

apoe_apgen <- dplyr::select(apoe_alt, rid, viscode, apgen1, apgen2)
apoe_apgen_sc <- dplyr::filter(apoe_apgen, viscode == "sc")
apoe_apgen_sc %>% count(!is.na(apgen2)) # TRUE 822
apoe_apgen_bl <- dplyr::filter(apoe_apgen, viscode == "bl")
apoe_apgen_bl %>% count(!is.na(apgen2)) # 0
apoe_apgen %>% count(!is.na(apgen2)) # TRUE 2556

# Join these two datasets, adsp & apgen
apoe_apgen_redef <- apoe_apgen %>%
  mutate(
    apoe4 = case_when(
      apgen1 == 4 & apgen2 == 4 ~ 2,
      (apgen1 == 4 & apgen2 != 4) | 
        (apgen1 != 4 & apgen2 == 4) ~ 1,
      apgen1 != 4 & apgen2 != 4 ~ 0,
      TRUE ~ NA_integer_
    ),
    apoe3 = case_when(
      apgen1 == 3 & apgen2 == 3 ~ 2,
      (apgen1 == 3 & apgen2 != 3) | 
        (apgen1 != 3 & apgen2 == 3) ~ 1,
      apgen1 != 3 & apgen2 != 3 ~ 0,
      TRUE ~ NA_integer_
    ),
    apoe2 = case_when(
      apgen1 == 2 & apgen2 == 2 ~ 2,
      (apgen1 == 2 & apgen2 != 2) | 
        (apgen1 != 2 & apgen2 == 2) ~ 1,
      apgen1 != 2 & apgen2 != 2 ~ 0,
      TRUE ~ NA_integer_
    )
  )

apoe_adsp <- dplyr::select(apoe_meta_name, rid, apoe4count, 
                           apoe3count, apoe2count)
apoe_adsp %>% count(!is.na(apoe4count)) # TRUE 7121, FALSE 3688
apoe_adsp_rename = dplyr::rename(apoe_adsp, apoe4 = apoe4count,
                                 apoe3 = apoe3count, apoe2 = apoe2count)
apoe_apgen_no_vis <- dplyr::select(apoe_apgen_redef, rid, apoe2, apoe3, 
                                   apoe4)
apoe_apgen_no_vis %>% count(!is.na(apoe4)) # TRUE 2556, FALSE 0
apoe_list <- list(apoe_adsp_rename, apoe_apgen_no_vis)
apoe_joined <- apoe_list %>% reduce(full_join, by=c('rid', 'apoe2'))

# don't use apoe_joined, or apoe_adsp; use 'apoe_apgen_redef'
str(apoe_apgen_redef) # rows 2556
apoe_apgen_redef[apoe_apgen_redef == -1] <- NA_real_
apoe_apgen_redef[apoe_apgen_redef == -4] <- NA_real_
str(apoe_apgen_redef) # rows 2556

str(apoe_joined) # 11816 rows
apoe_final <- apoe_apgen_redef %>% dplyr::distinct()
str(apoe_final) # 2556 rows
apoe_final %>% count(!is.na(apoe2)) # TRUE 2556, FALSE 0

############### APOE computed #####################
str(imputed_wo_chol_lj) # rows 2430
summary(data_final) # adas11_bl na: 9, m_pacc na: 2
# Summary so far: 
# We have apoe (apoe_final): Contains rid, apoe2, apoe3, apoe4
# We have caide vars imputed (imputed_wo_chol_lj): 
#         Contains rid, age, edu (pteducat), sex (sex_factor), 
#                  bmi, sysbp (vsbpsys)
# We have non-caide vars non-imputed (data_final): Contains:
#         rid, race, dx_bl, adas11_bl, m_pac_cdigit_bl
# We have chol (chol_final): Contains: rid, total_c

# Need to make two datasets:
# 1. dataset_w_chol for caide_calc; then drop people wo chol (~1500 people)
# 2. dataset_wo_chol for caide_calc; don't drop anybody (~2200 people)

# Join datasets:
# Dataset 1:
impute_list1 <- list(apoe_final, imputed_wo_chol, data_final, chol_final)
joined_w_chol <- impute_list1 %>% reduce(full_join, by='rid')

# Dataset 2:
impute_list2 <- list(apoe_final, imputed_wo_chol_lj, data_final)
joined_wo_chol <- impute_list2 %>% reduce(full_join, by='rid')

######################################################
# left_join to make the complete dataset
######################################################
summary(imputed_wo_chol_lj)
summary(data_final)
apoe_final %>% count(!(is.na(apoe2))) # na 0

joined_wo_chol0 <- imputed_wo_chol_lj %>% left_join(data_final, 
                                              by=c('rid'))
joined_wo_chol1 <- joined_wo_chol0 %>% left_join(apoe_final, 
                                                   by=c('rid'))

str(joined_wo_chol1) #2430
joined_wo_chol4 <- dplyr::select(joined_wo_chol1, rid, age, pteducat, 
                                sex_factor, bmi, vsbpsys, race, 
                                dx_bl, apoe4, apoe3, apoe2)
str(joined_wo_chol4) #2430
joined_wo_chol <- na.omit(joined_wo_chol4)
str(joined_wo_chol) # rows 2174 removed adas & mpacc

joined_wo_chol %>% count(apoe2) # na 0
# previous;y, checked: race na 61, dx_bl 61, adas11_bl 70, m_pac_cdigit_bl 63

# This is: Left_join with imputed_wo_chol_lj


########################################################
# Dataset 1 calculations: dataset_w_chol for caide_calc; 
#                         then drop people wo chol (~1500 people)
########################################################

# Compute apoe_allele & hypchol:
data_convert <- joined_w_chol %>%
  mutate(
    hypchol = case_when(
      total_c > 6.18 ~ 1,
      total_c <= 6.18 ~ 0,
      TRUE ~ NA_real_
    ),
    apoe_allele = case_when(
      apoe4 == 1 & apoe2 == 1 ~ 24,
      apoe4 == 1 & apoe3 == 1 ~ 34,
      apoe4 == 0 & apoe2 == 2 ~ 22,
      apoe4 == 0 & apoe3 == 2 ~ 33,
      apoe4 == 2 & apoe2 == 0 ~ 44,
      apoe4 == 2 & apoe3 == 0 ~ 44,
      apoe2 == 1 & apoe3 == 1 ~ 23,
      apoe2 == 2 & apoe3 == 0 ~ 22,
      apoe3 == 2 & apoe2 == 0 ~ 33,
      TRUE ~ NA_real_
    )
  )

# Relabel the imputed variables:
joined_w_chol_i = dplyr::rename(data_convert, i_age = age, i_education =
                                  pteducat, i_bmi = bmi, i_systbp = vsbpsys)

# Compute caide with chol:
# 1. caide calculation
caide_calc <- joined_w_chol_i %>%
  mutate(
    caide_age = case_when(
      i_age < 46 ~ 0,
      i_age >= 46 & i_age <= 53 ~ 3,
      i_age > 53 ~ 4,
      TRUE ~ NA_real_
    ),
    caide_educ = case_when(
      i_education < 7 ~ 3,
      i_education >= 7 & i_education <= 9 ~ 2,
      i_education > 9 ~ 0,
      TRUE ~ NA_real_
    ),
    caide_sex = case_when(
      sex_factor == "Male" ~ 1,
      sex_factor == "Female" ~ 0,
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
    caide_chol = case_when(
      hypchol == 0 ~ 0, 
      hypchol == 1 ~ 2, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    caide = sum(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
                caide_chol, na.rm = F),
    caide_missing = sum(is.na(i_age), is.na(i_education), is.na(sex_factor), 
                        is.na(i_bmi), is.na(i_systbp), is.na(hypchol))
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
      apoe_allele == 22 ~ 'e2+',
      apoe_allele == 23 ~ 'e2+',
      apoe_allele == 24 ~ 'e4+',
      apoe_allele == 33 ~ 'e3/e3',
      apoe_allele == 34 ~ 'e4+',
      apoe_allele == 44 ~ 'e4+',
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
         caide_chol), as_factor
  ) %>%
  dplyr::select(rid, starts_with('caide'), z_caide, caide_cat, apoe, race,
                caide_apoe, i_age, i_education, sex_factor, i_bmi, i_systbp, 
                hypchol, caide_chol, dx_bl, adas11_bl, 
                m_pac_cdigit_bl)

write_csv(x = caide_calc, "caide_adni_chol.csv") 

####################################################
# Date: Wed, Sept 6, 2023, 9:13am
# finish: need mcaide calculation
# Need dataset 2 & mcaide calc: Work on dataset 2 first
####################################################

caide_calc <- read_csv("caide_adni_chol.csv") %>% 
  janitor::clean_names()

str(caide_calc) # rows = 3599; na.omit & see how many rows survive
caide_calc_na_omit_all <- as.data.frame(caide_calc)
caide_na_omit <- na.omit(caide_calc_na_omit_all)
str(caide_na_omit) # rows = 1691

summary(caide_calc_na_omit_all)

########################################################
# Dataset2 calculation:
# caide calculation without chol; do not drop people (~2200)
########################################################

# Compute caide without chol:
# 1. caide calculation
# Dataset2: joined_wo_chol

# Compute apoe_allele:
data_convert_c <- joined_wo_chol %>%
  mutate(
    apoe_allele = case_when(
      apoe4 == 1 & apoe2 == 1 ~ 24,
      apoe4 == 1 & apoe3 == 1 ~ 34,
      apoe4 == 0 & apoe2 == 2 ~ 22,
      apoe4 == 0 & apoe3 == 2 ~ 33,
      apoe4 == 2 & apoe2 == 0 ~ 44,
      apoe4 == 2 & apoe3 == 0 ~ 44,
      apoe2 == 1 & apoe3 == 1 ~ 23,
      apoe2 == 2 & apoe3 == 0 ~ 22,
      apoe3 == 2 & apoe2 == 0 ~ 33,
      TRUE ~ NA_real_
    )
  )

# Relabel the imputed variables:
joined_wo_chol_i = dplyr::rename(data_convert_c, i_age = age, i_education =
                                  pteducat, i_bmi = bmi, i_systbp = vsbpsys)

caide_calc_nc <- joined_wo_chol_i %>%
  mutate(
    caide_age = case_when(
      i_age < 46 ~ 0,
      i_age >= 46 & i_age <= 53 ~ 3,
      i_age > 53 ~ 4,
      TRUE ~ NA_real_
    ),
    caide_educ = case_when(
      i_education < 7 ~ 3,
      i_education >= 7 & i_education <= 9 ~ 2,
      i_education > 9 ~ 0,
      TRUE ~ NA_real_
    ),
    caide_sex = case_when(
      sex_factor == "Male" ~ 1,
      sex_factor == "Female" ~ 0,
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
    )
  ) %>%
  rowwise() %>%
  mutate(
    caide = sum(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
                 na.rm = F),
    caide_missing = sum(is.na(i_age), is.na(i_education), is.na(sex_factor), 
                        is.na(i_bmi), is.na(i_systbp))
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
      apoe_allele == 22 ~ 'e2+',
      apoe_allele == 23 ~ 'e2+',
      apoe_allele == 24 ~ 'e4+',
      apoe_allele == 33 ~ 'e3/e3',
      apoe_allele == 34 ~ 'e4+',
      apoe_allele == 44 ~ 'e4+',
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
    vars(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp 
         ), as_factor
  ) %>%
  dplyr::select(rid, starts_with('caide'), z_caide, caide_cat, apoe, race,
                caide_apoe, i_age, i_education, sex_factor, i_bmi, i_systbp, 
                dx_bl)

write_csv(x = caide_calc_nc, "caide_adni_nc.csv") 

# 2. Compute mcaide now (without chol):

data_convert_m <- joined_wo_chol %>%
  mutate(
    apoe_allele = case_when(
      apoe4 == 1 & apoe2 == 1 ~ 24,
      apoe4 == 1 & apoe3 == 1 ~ 34,
      apoe4 == 0 & apoe2 == 2 ~ 22,
      apoe4 == 0 & apoe3 == 2 ~ 33,
      apoe4 == 2 & apoe2 == 0 ~ 44,
      apoe4 == 2 & apoe3 == 0 ~ 44,
      apoe2 == 1 & apoe3 == 1 ~ 23,
      apoe2 == 2 & apoe3 == 0 ~ 22,
      apoe3 == 2 & apoe2 == 0 ~ 33,
      TRUE ~ NA_real_
    )
  )

# Relabel the imputed variables:
joined_wo_chol_i = dplyr::rename(data_convert_m, i_age = age, i_education =
                                  pteducat, i_bmi = bmi, i_systbp = vsbpsys)

mcaide_calc_nc <- joined_wo_chol_i %>%
  mutate(
    mcaide_age = case_when(
      i_age < 64 ~ 0,
      i_age >= 64 & i_age <= 73 ~ 1,
      i_age > 73 ~ 2,
      TRUE ~ NA_real_
    ),
    mcaide_educ = case_when(
      i_education < 12 ~ 2,
      i_education >= 12 & i_education <= 16 ~ 1,
      i_education > 16 ~ 0,
      TRUE ~ NA_real_
    ),
    mcaide_sex = case_when(
      sex_factor == "Male" ~ 1,
      sex_factor == "Female" ~ 0,
      TRUE ~ NA_real_
    ),
    mcaide_obesity = case_when(
      i_bmi < 30 ~ 0,
      i_bmi >= 30 ~ 2,
      TRUE ~ NA_real_
    ),
    mcaide_sbp = case_when(
      i_systbp < 140 ~ 0,
      i_systbp >= 140 ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    mcaide = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, 
                 mcaide_sbp, na.rm = F),
    mcaide_missing = sum(is.na(i_age), is.na(i_education), is.na(sex_factor), 
                         is.na(i_bmi), is.na(i_systbp))
  ) %>%
  ungroup() %>%
  mutate(
    z_mcaide = scale(mcaide)[,1],
    mcaide_cat = case_when(
      between(mcaide, 0, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T))) 
      ~ 'Low',
      between(mcaide, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T)), 
              (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T))) ~ 'Mid',
      between(mcaide, (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    apoe = case_when(
      apoe_allele == 22 ~ 'e2+',
      apoe_allele == 23 ~ 'e2+',
      apoe_allele == 24 ~ 'e4+',
      apoe_allele == 33 ~ 'e3/e3',
      apoe_allele == 34 ~ 'e4+',
      apoe_allele == 44 ~ 'e4+',
      TRUE ~ 'NA'
    )
  ) %>%
  mutate(
    mcaide_cat = fct_relevel(mcaide_cat, 'Low', "Mid",  'High'),
    mcaide_apoe = glue("{mcaide_cat}, {apoe}"),
    mcaide_apoe = fct_relevel(mcaide_apoe, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                              "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                              "High, e4+")
  ) %>%
  mutate_at(
    vars(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp), as_factor
  ) %>%
  dplyr::select(rid, starts_with('mcaide'), z_mcaide, mcaide_cat, apoe, race,
                mcaide_apoe, i_age, i_education, sex_factor, i_bmi, i_systbp, 
                dx_bl)

write_csv(x = mcaide_calc_nc, "mcaide_adni_nc.csv")

# check # of rows:
str(mcaide_calc_nc) # 2174 rows
str(caide_calc_nc) # 2174 rows

no_na_caide <- na.omit(caide_calc_nc)
no_na_mcaide <- na.omit(mcaide_calc_nc)
str(no_na_caide) # rows = 2408
str(no_na_mcaide) # rows = 2408

##################################################################
# Make tables for caide and mcaide (without chol)
# Tables are stratified by RACE
# 
# Dataset #2:
# caide_calc_nc # 2174 rows
# mcaide_calc_nc # 2174 rows

##################################################################

# 1. caide table

caide_1 <- caide_calc_nc %>%
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
    )
  )

caide_1 %>%
  dplyr::select(dx_bl, race,
                Gender, Age, Education, Obesity, Hypertension, 
                caide, apoe, caide_apoe 
  ) %>%
  tbl_summary(
    by = race,
    type = list(caide ~ "continuous", Gender ~ "categorical", 
                Obesity ~ "categorical", Hypertension ~ "categorical", 
                apoe ~ "categorical", race ~ "categorical",
                dx_bl ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(Gender ~ "Gender", Age ~ "Age", Education ~ "Education", 
                 race ~ "Race", dx_bl ~ "Diagnosis",
                 Obesity ~ "Obesity", Hypertension ~ "Hypertension", 
                  caide ~ "CAIDE", 
                 apoe ~ "APOE", caide_apoe ~ "CAIDE x APOE"),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**RACE**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "table_caide_adni_wo_chol.html"
  )

webshot("table_caide_adni_wo_chol.html", "table_caide_adni_wo_chol.png")

# 2. mcaide table

mcaide_1 <- mcaide_calc_nc %>%
  mutate(
    Gender = case_when(
      (mcaide_sex == 1) ~ "Male",
      (mcaide_sex == 0) ~ "Female",
      is.na(mcaide_sex) ~ NA
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
    filename = "table_mcaide_adni_wo_chol.html"
  )

webshot("table_mcaide_adni_wo_chol.html", "table_mcaide_adni_wo_chol.png")

##########################################################
# Regression and Odds Ratio
# m/caide dataset without cholesterol in the calculation
##########################################################

adni_caide_binary <- caide_calc_nc %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary <- mcaide_calc_nc %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ m/caide + apoe
adni_caide_binary$apoe <- factor(adni_caide_binary$apoe, ordered = FALSE )
adni_caide_binary$apoe = relevel(adni_caide_binary$apoe, ref = "e3/e3")
adni_mcaide_binary$apoe <- factor(adni_mcaide_binary$apoe, ordered = FALSE )
adni_mcaide_binary$apoe = relevel(adni_mcaide_binary$apoe, ref = "e3/e3")

# 1. logistic model: diagnosis ~ caide + apoe (without cholesterol)
model_caide_logistic <- glm(diag ~ caide + apoe, 
                                 data=adni_caide_binary, family = binomial)

# 2. logistic model: diagnosis ~ mcaide + apoe (without cholesterol)
model_mcaide_logistic <- glm(diag ~ mcaide + apoe, 
                                  data=adni_mcaide_binary, family = binomial)

# 3 logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education

model_caide_logistic_all <- glm(diag ~ i_bmi + i_systbp + 
                                        i_age + i_education + sex_factor + 
                                  apoe, 
                                      data=adni_caide_binary, 
                                family = binomial)

# 4 logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education +
#                               sex_factor + apoe
model_mcaide_logistic_all <- glm(diag ~ i_bmi + i_systbp + 
                                         i_age + i_education + sex_factor + apoe, 
                                       data=adni_mcaide_binary, 
                                       family = binomial)

##########################################################
# Glance
##########################################################

add_glance_source_note(tbl_regression(model_mcaide_logistic_all)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table. Logistic Regression for mCAIDE, APOE") %>%
  gt::gtsave(
    filename = "model_mcaide_logistic_all.html"
  )

webshot("model_mcaide_logistic_all.html", "model_mcaide_logistic_all.png")

##########################################################
# Compute Odds Ratio:
##########################################################

# results give estimate and std.error:
results <- summary(model_caide_logistic)$coefficients
df_apoe_c <- as.data.frame(results)
estimate <- df_apoe_c$Estimate[-1]
std_err <- df_apoe_c$`Std. Error`[-1]

apoe <- c("CAIDE", "e2+", "e4+")
data <- tibble(apoe, estimate, std_err)
colnames(data) <- c("apoe", "estimate", "std_err")

data_or <- data %>% mutate(
  lci = estimate - (std_err * 1.96),
  hci = estimate + (std_err * 1.96), 
  or = exp(estimate),
  or_lci = exp(lci),
  or_hci = exp(hci),
  colorx = case_when(
    (or_lci > 1) | (or_hci < 1) ~ "red"
  )
)

ggplot(data_or, aes(x = or, y = factor(apoe, levels = c('CAIDE', 'e2+', 
                                                        'e4+'),
), color = colorx)
) +
  #ggplot(data_or, aes(x = or, y = factor(apoe, levels = c("Low, e2+", 
  #          "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+", "High, e2+", 
  #          "High, e3/e3", "High, e4+")
  #          ))) +
  labs(title = "CAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci, xmax = or_hci), linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.01)) + 
  # scale_color_manual(values = c( "#E41A1C", "#377EB8", "#4DAF4A")) +
  theme_bw() + 
  theme(
    text = element_text(size = 15),
    # axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5)
  )

########################################################
# Regression stratified by race
########################################################
write_csv(caide_calc_nc, "caide_calc_nc.csv")
write_csv(mcaide_calc_nc, "mcaide_calc_nc.csv")

caide_calc_nc <- read_csv("caide_calc_nc.csv") %>% 
  janitor::clean_names()

mcaide_calc_nc <- read_csv("mcaide_calc_nc.csv") %>% 
  janitor::clean_names()

race_asian_caide <-  dplyr::filter(caide_calc_nc, 
                                 race == "Asian")
race_asian_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                                  race == "Asian")

adni_caide_binary <- race_asian_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary <- race_asian_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ m/caide + apoe
adni_caide_binary$apoe <- factor(adni_caide_binary$apoe, ordered = FALSE )
adni_caide_binary$apoe = relevel(adni_caide_binary$apoe, ref = "e3/e3")
adni_mcaide_binary$apoe <- factor(adni_mcaide_binary$apoe, ordered = FALSE )
adni_mcaide_binary$apoe = relevel(adni_mcaide_binary$apoe, ref = "e3/e3")

# regression models: FOR asian race:
# 1. logistic model: diagnosis ~ caide + apoe (without cholesterol)
model_caide_logistic_asian <- glm(diag ~ caide + apoe, 
                            data=adni_caide_binary, family = binomial)

# 2. logistic model: diagnosis ~ mcaide + apoe (without cholesterol)
model_mcaide_logistic_asian <- glm(diag ~ mcaide + apoe, 
                             data=adni_mcaide_binary, family = binomial)

# 3 logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education
              
model_caide_logistic_all_asian <- glm(diag ~ i_bmi + i_systbp + 
                              i_age + i_education + sex_factor + apoe, 
                              data=adni_caide_binary, family = binomial)

# 4 logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education +
#                               sex_factor + apoe
model_mcaide_logistic_all_asian <- glm(diag ~ i_bmi + i_systbp + 
                                    i_age + i_education + sex_factor + apoe, 
                                    data=adni_mcaide_binary, 
                                    family = binomial)

#############################################
# Glance
#############################################

# 1. model for asian & caide
add_glance_source_note(tbl_regression(model_caide_logistic_asian)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table. Logistic Regression for 
                 CAIDE (Asian)") %>%
  gt::gtsave(
    filename = "model_caide_logistic_asian.html"
  )

webshot("model_caide_logistic_asian.html", "model_caide_logistic_asian.png")

# 2. model for asian & mcaide
add_glance_source_note(tbl_regression(model_mcaide_logistic_asian)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table. Logistic Regression for 
                 mCAIDE (Asian)") %>%
  gt::gtsave(
    filename = "model_mcaide_logistic_asian.html"
  )

webshot("model_mcaide_logistic_asian.html", "model_mcaide_logistic_asian.png")

# 3. model for asian & caide & all continuous vars
add_glance_source_note(tbl_regression(model_caide_logistic_all_asian)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table. Logistic Regression for 
                 CAIDE (Asian)") %>%
  gt::gtsave(
    filename = "model_caide_logistic_all_asian.html"
  )

webshot("model_caide_logistic_all_asian.html", 
        "model_caide_logistic_all_asian.png")

# 4. model for asian & mcaide & all continuous vars
add_glance_source_note(tbl_regression(model_mcaide_logistic_all_asian)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table. Logistic Regression for 
          mCAIDE (Asian)") %>%
  gt::gtsave(
    filename = "model_mcaide_logistic_all_asian.html"
  )

webshot("model_mcaide_logistic_all_asian.html", 
        "model_mcaide_logistic_all_asian.png")

############################################################
# mCAIDE: Odds Ratio Figures for Slide 28: Stratified by Race
############################################################

mcaide_calc_nc <- read_csv("mcaide_calc_nc.csv") %>% 
  janitor::clean_names()

nhw_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                                   race == "Non-Hispanic White")
asian_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                                  race == "Asian")
black_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                               race == "Black")
hispanic_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                               race == "Hispanic")

adni_mcaide_binary_nhw <- nhw_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_asian <- asian_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_black <- black_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_his <- hispanic_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ mcaide + apoe
adni_mcaide_binary_nhw$apoe <- factor(adni_mcaide_binary_nhw$apoe, 
                                      ordered = FALSE )
adni_mcaide_binary_nhw$apoe = relevel(adni_mcaide_binary_nhw$apoe, 
                                      ref = "e3/e3")
adni_mcaide_binary_asian$apoe <- factor(adni_mcaide_binary_asian$apoe, 
                                     ordered = FALSE )
adni_mcaide_binary_asian$apoe = relevel(adni_mcaide_binary_asian$apoe, 
                                     ref = "e3/e3")
adni_mcaide_binary_black$apoe <- factor(adni_mcaide_binary_black$apoe, 
                                        ordered = FALSE )
adni_mcaide_binary_black$apoe = relevel(adni_mcaide_binary_black$apoe, 
                                        ref = "e3/e3")
adni_mcaide_binary_his$apoe <- factor(adni_mcaide_binary_his$apoe, 
                                        ordered = FALSE )
adni_mcaide_binary_his$apoe = relevel(adni_mcaide_binary_his$apoe, 
                                        ref = "e3/e3")

# 1. logistic model: diagnosis ~ caide + apoe (without cholesterol)
model_mcaide_logistic_nhw <- glm(diag ~ mcaide + apoe, 
                            data=adni_mcaide_binary_nhw, 
                            family = binomial)
model_mcaide_logistic_asian <- glm(diag ~ mcaide + apoe, 
                                 data=adni_mcaide_binary_asian, 
                                 family = binomial)
model_mcaide_logistic_black <- glm(diag ~ mcaide + apoe, 
                                   data=adni_mcaide_binary_black, 
                                   family = binomial)
model_mcaide_logistic_his <- glm(diag ~ mcaide + apoe, 
                                   data=adni_mcaide_binary_his, 
                                   family = binomial)


# results give estimate and std.error for e2+:
results_nhw <- summary(model_mcaide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_e2 <- df_apoe_nhw$Estimate[3]
std_err_nhw_e2 <- df_apoe_nhw$`Std. Error`[3]

results_asian <- summary(model_mcaide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_e2 <- df_apoe_asian$Estimate[3]
std_err_asian_e2 <- df_apoe_asian$`Std. Error`[3]

results_black <- summary(model_mcaide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_e2 <- df_apoe_black$Estimate[3]
std_err_black_e2 <- df_apoe_black$`Std. Error`[3]

results_hispanic <- summary(model_mcaide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_e2 <- df_apoe_his$Estimate[3]
std_err_his_e2 <- df_apoe_his$`Std. Error`[3]

estimate_e2 <- c(estimate_nhw_e2, estimate_asian_e2, estimate_black_e2, 
              estimate_his_e2)
std_err_e2 <- c(std_err_nhw_e2, std_err_asian_e2, std_err_black_e2, 
             std_err_his_e2)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
data_e2 <- tibble(apoe, estimate_e2, std_err_e2)
colnames(data_e2) <- c("apoe", "estimate_e2", "std_err_e2")

data_or_e2 <- data_e2 %>% mutate(
  lci_e2 = estimate_e2 - (std_err_e2 * 1.96),
  hci_e2 = estimate_e2 + (std_err_e2 * 1.96), 
  or_e2 = exp(estimate_e2),
  or_lci_e2 = exp(lci_e2),
  or_hci_e2 = exp(hci_e2),
  colorx = case_when(
    (or_lci_e2 > 1) | (or_hci_e2 < 1) ~ "red"
  )
)

# odds ratio/coefficient figure"
p1 <- ggplot(data_or_e2, aes(x = or_e2, y = factor(apoe, levels = c('NHW', 'Asian', 
                                                        'Black', 'Hispanic')), 
                    color = colorx)
) +
  #ggplot(data_or, aes(x = or, y = factor(apoe, levels = c("Low, e2+", 
  #          "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+", "High, e2+", 
  #          "High, e3/e3", "High, e4+")
  #          ))) +
  labs(title = "e2+") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci_e2, xmax = or_hci_e2), linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous(" ", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  # scale_color_manual(values = c( "#E41A1C", "#377EB8", "#4DAF4A")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0,-5,0,-2), "cm")
  )

# results give estimate and std.error for e4+:
results_nhw <- summary(model_mcaide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_e4 <- df_apoe_nhw$Estimate[4]
std_err_nhw_e4 <- df_apoe_nhw$`Std. Error`[4]

results_asian <- summary(model_mcaide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_e4 <- df_apoe_asian$Estimate[4]
std_err_asian_e4 <- df_apoe_asian$`Std. Error`[4]

results_black <- summary(model_mcaide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_e4 <- df_apoe_black$Estimate[4]
std_err_black_e4 <- df_apoe_black$`Std. Error`[4]

results_hispanic <- summary(model_mcaide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_e4 <- df_apoe_his$Estimate[4]
std_err_his_e4 <- df_apoe_his$`Std. Error`[4]

estimate_e4 <- c(estimate_nhw_e4, estimate_asian_e4, estimate_black_e4, 
              estimate_his_e4)
std_err_e4 <- c(std_err_nhw_e4, std_err_asian_e4, std_err_black_e4, 
             std_err_his_e4)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
data_e4 <- tibble(apoe, estimate_e4, std_err_e4)
colnames(data_e4) <- c("apoe", "estimate_e4", "std_err_e4")

data_or_e4 <- data_e4 %>% mutate(
  lci_e4 = estimate_e4 - (std_err_e4 * 1.96),
  hci_e4 = estimate_e4 + (std_err_e4 * 1.96), 
  or_e4 = exp(estimate_e4),
  or_lci_e4 = exp(lci_e4),
  or_hci_e4 = exp(hci_e4),
  colorx = case_when(
    (or_lci_e4 > 1) | (or_hci_e4 < 1) ~ "red"
  )
)

# odds ratio/coefficient figure"
p2 <- ggplot(data_or_e4, aes(x = or_e4, y = factor(apoe, levels = c('NHW', 'Asian', 
                                                        'Black', 'Hispanic')), 
                    color = colorx)
) +
  #ggplot(data_or, aes(x = or, y = factor(apoe, levels = c("Low, e2+", 
  #          "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+", "High, e2+", 
  #          "High, e3/e3", "High, e4+")
  #          ))) +
  labs(title = "e4+") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci_e4, xmax = or_hci_e4), linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_color_manual(values = c( "#E41A1C", "#377EB8", "#4DAF4A")) +
  theme_bw() + 
  theme(
    # axis.title=element_text(size=12,face="bold"),
    axis.title=element_text(size=12),
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# results give estimate and std.error for "mcaide":
results_nhw <- summary(model_mcaide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_mcaide <- df_apoe_nhw$Estimate[2]
std_err_nhw_mcaide <- df_apoe_nhw$`Std. Error`[2]

results_asian <- summary(model_mcaide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_mcaide <- df_apoe_asian$Estimate[2]
std_err_asian_mcaide <- df_apoe_asian$`Std. Error`[2]

results_black <- summary(model_mcaide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_mcaide <- df_apoe_black$Estimate[2]
std_err_black_mcaide <- df_apoe_black$`Std. Error`[2]

results_hispanic <- summary(model_mcaide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_mcaide <- df_apoe_his$Estimate[2]
std_err_his_mcaide <- df_apoe_his$`Std. Error`[2]

estimate_mcaide <- c(estimate_nhw_mcaide, estimate_asian_mcaide, 
              estimate_black_mcaide, 
              estimate_his_mcaide)
std_err_mcaide <- c(std_err_nhw_mcaide, std_err_asian_mcaide, std_err_black_mcaide, 
             std_err_his_mcaide)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
data_mcaide <- tibble(apoe, estimate_mcaide, std_err_mcaide)
colnames(data_mcaide) <- c("apoe", "estimate_mcaide", "std_err_mcaide")

data_or_mcaide <- data_mcaide %>% mutate(
  lci_mcaide = estimate_mcaide - (std_err_mcaide * 1.96),
  hci_mcaide = estimate_mcaide + (std_err_mcaide * 1.96), 
  or_mcaide = exp(estimate_mcaide),
  or_lci_mcaide = exp(lci_mcaide),
  or_hci_mcaide = exp(hci_mcaide),
  colorx = case_when(
    (or_lci_mcaide > 1) | (or_hci_mcaide < 1) ~ "red"
  )
) 

# odds ratio/coefficient figure"
p3 <- ggplot(data_or_mcaide, aes(x = or_mcaide, y = factor(apoe, 
                                 levels = c('NHW', 'Asian', 
                                 'Black', 'Hispanic')),
                                 color = colorx)
) +
  #ggplot(data_or, aes(x = or, y = factor(apoe, levels = c("Low, e2+", 
  #          "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+", "High, e2+", 
  #          "High, e3/e3", "High, e4+")
  #          ))) +
  labs(title = "mCAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci_mcaide, xmax = or_hci_mcaide), 
                 linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous(" ", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  # scale_color_manual(values = c( "#E41A1C", "#377EB8", "#4DAF4A")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0,-1.5,0,-2), "cm")
  )

# fix color error:
data_or_mcaide$colorx[is.na(data_or_mcaide$colorx)] <- "black"
data_or_e2$colorx[is.na(data_or_e2$colorx)] <- "black"
data_or_e4$colorx[is.na(data_or_e4$colorx)] <- "black"


p1 + p2 + p3 + plot_layout(nrow = 1, ncol = 3, heights = c(0.5, 1))

############################################################
# CAIDE: Odds Ratio Figures for Slide 28: Facets
############################################################
# fix color error:
data_or_mcaide$colorx[is.na(data_or_mcaide$colorx)] <- "black"
data_or_e2$colorx[is.na(data_or_e2$colorx)] <- "black"
data_or_e4$colorx[is.na(data_or_e4$colorx)] <- "black"


var <- c("mCAIDE", "mCAIDE", "mCAIDE", "mCAIDE")
mcaide_facets1 = tibble(data_or_mcaide, var)

var <- c("e2+", "e2+", "e2+", "e2+")
e2_facets1 = tibble(data_or_e2, var)

var <- c("e4+", "e4+", "e4+", "e4+")
e4_facets1 = tibble(data_or_e4, var)


# rename and concatenate:
mcaide_facets <- dplyr::rename(mcaide_facets1, estimate = estimate_mcaide,
                               std_err = std_err_mcaide, lci = lci_mcaide,
                               hci = hci_mcaide, or = or_mcaide, race = apoe,
                               or_lci = or_lci_mcaide, or_hci = or_hci_mcaide)
e2_facets <- dplyr::rename(e2_facets1, estimate = estimate_e2,
                           std_err = std_err_e2, lci = lci_e2,
                           hci = hci_e2, or = or_e2, or_lci = or_lci_e2,
                           or_hci = or_hci_e2, race = apoe)
e4_facets <- dplyr::rename(e4_facets1, estimate = estimate_e4,
                           std_err = std_err_e4, lci = lci_e4,
                           hci = hci_e4, or = or_e4, or_lci = or_lci_e4,
                           or_hci = or_hci_e4, race = apoe)

data_facets <- rbind(e2_facets, e4_facets, mcaide_facets)

ggplot(data_facets, aes(x=or, y = factor(race, levels = c('NHW', 'Asian', 
                        'Black', 'Hispanic')), group=var, color = colorx)) +
  facet_grid(. ~ var, scales='free_x') +
  #labs(title = "mCAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci, xmax = or_hci), 
                 linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c( "black", "red")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5)
  )
        
############################################################
# CAIDE: Odds Ratio Figures for Slide 28: Stratified by Race
############################################################

caide_calc_nc <- read_csv("caide_calc_nc.csv") %>% 
  janitor::clean_names()

nhw_caide <-  dplyr::filter(caide_calc_nc, 
                             race == "Non-Hispanic White")
asian_caide <-  dplyr::filter(caide_calc_nc, 
                               race == "Asian")
black_caide <-  dplyr::filter(caide_calc_nc, 
                               race == "Black")
hispanic_caide <-  dplyr::filter(caide_calc_nc, 
                                  race == "Hispanic")

adni_caide_binary_nhw <- nhw_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_asian <- asian_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_black <- black_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_his <- hispanic_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ mcaide + apoe
adni_caide_binary_nhw$apoe <- factor(adni_caide_binary_nhw$apoe, 
                                      ordered = FALSE )
adni_caide_binary_nhw$apoe = relevel(adni_caide_binary_nhw$apoe, 
                                      ref = "e3/e3")
adni_caide_binary_asian$apoe <- factor(adni_caide_binary_asian$apoe, 
                                        ordered = FALSE )
adni_caide_binary_asian$apoe = relevel(adni_caide_binary_asian$apoe, 
                                        ref = "e3/e3")
adni_caide_binary_black$apoe <- factor(adni_caide_binary_black$apoe, 
                                        ordered = FALSE )
adni_caide_binary_black$apoe = relevel(adni_caide_binary_black$apoe, 
                                        ref = "e3/e3")
adni_caide_binary_his$apoe <- factor(adni_caide_binary_his$apoe, 
                                      ordered = FALSE )
adni_caide_binary_his$apoe = relevel(adni_caide_binary_his$apoe, 
                                      ref = "e3/e3")

# 1. logistic model: diagnosis ~ caide + apoe (without cholesterol)
model_caide_logistic_nhw <- glm(diag ~ caide + apoe, 
                                 data=adni_caide_binary_nhw, 
                                 family = binomial)
model_caide_logistic_asian <- glm(diag ~ caide + apoe, 
                                   data=adni_caide_binary_asian, 
                                   family = binomial)
model_caide_logistic_black <- glm(diag ~ caide + apoe, 
                                   data=adni_caide_binary_black, 
                                   family = binomial)
model_caide_logistic_his <- glm(diag ~ caide + apoe, 
                                 data=adni_caide_binary_his, 
                                 family = binomial)


# results give estimate and std.error for e2+:
results_nhw <- summary(model_caide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_e2 <- df_apoe_nhw$Estimate[3]
std_err_nhw_e2 <- df_apoe_nhw$`Std. Error`[3]

results_asian <- summary(model_caide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_e2 <- df_apoe_asian$Estimate[3]
std_err_asian_e2 <- df_apoe_asian$`Std. Error`[3]

results_black <- summary(model_caide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_e2 <- df_apoe_black$Estimate[3]
std_err_black_e2 <- df_apoe_black$`Std. Error`[3]

results_hispanic <- summary(model_caide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_e2 <- df_apoe_his$Estimate[3]
std_err_his_e2 <- df_apoe_his$`Std. Error`[3]

estimate_e2 <- c(estimate_nhw_e2, estimate_asian_e2, estimate_black_e2, 
                 estimate_his_e2)
std_err_e2 <- c(std_err_nhw_e2, std_err_asian_e2, std_err_black_e2, 
                std_err_his_e2)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
data_e2 <- tibble(apoe, estimate_e2, std_err_e2)
colnames(data_e2) <- c("apoe", "estimate_e2", "std_err_e2")

data_or_e2 <- data_e2 %>% mutate(
  lci_e2 = estimate_e2 - (std_err_e2 * 1.96),
  hci_e2 = estimate_e2 + (std_err_e2 * 1.96), 
  or_e2 = exp(estimate_e2),
  or_lci_e2 = exp(lci_e2),
  or_hci_e2 = exp(hci_e2),
  colorx = case_when(
    (or_lci_e2 > 1) | (or_hci_e2 < 1) ~ "red"
  )
)

# results give estimate and std.error for e4+:
results_nhw <- summary(model_caide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_e4 <- df_apoe_nhw$Estimate[4]
std_err_nhw_e4 <- df_apoe_nhw$`Std. Error`[4]

results_asian <- summary(model_caide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_e4 <- df_apoe_asian$Estimate[4]
std_err_asian_e4 <- df_apoe_asian$`Std. Error`[4]

results_black <- summary(model_caide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_e4 <- df_apoe_black$Estimate[4]
std_err_black_e4 <- df_apoe_black$`Std. Error`[4]

results_hispanic <- summary(model_caide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_e4 <- df_apoe_his$Estimate[4]
std_err_his_e4 <- df_apoe_his$`Std. Error`[4]

estimate_e4 <- c(estimate_nhw_e4, estimate_asian_e4, estimate_black_e4, 
                 estimate_his_e4)
std_err_e4 <- c(std_err_nhw_e4, std_err_asian_e4, std_err_black_e4, 
                std_err_his_e4)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
data_e4 <- tibble(apoe, estimate_e4, std_err_e4)
colnames(data_e4) <- c("apoe", "estimate_e4", "std_err_e4")

data_or_e4 <- data_e4 %>% mutate(
  lci_e4 = estimate_e4 - (std_err_e4 * 1.96),
  hci_e4 = estimate_e4 + (std_err_e4 * 1.96), 
  or_e4 = exp(estimate_e4),
  or_lci_e4 = exp(lci_e4),
  or_hci_e4 = exp(hci_e4),
  colorx = case_when(
    (or_lci_e4 > 1) | (or_hci_e4 < 1) ~ "red"
  )
)

# results give estimate and std.error for "caide":
results_nhw <- summary(model_caide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_caide <- df_apoe_nhw$Estimate[2]
std_err_nhw_caide <- df_apoe_nhw$`Std. Error`[2]

results_asian <- summary(model_caide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_caide <- df_apoe_asian$Estimate[2]
std_err_asian_caide <- df_apoe_asian$`Std. Error`[2]

results_black <- summary(model_caide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_caide <- df_apoe_black$Estimate[2]
std_err_black_caide <- df_apoe_black$`Std. Error`[2]

results_hispanic <- summary(model_caide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_hispanic)
estimate_his_caide <- df_apoe_his$Estimate[2]
std_err_his_caide <- df_apoe_his$`Std. Error`[2]

estimate_caide <- c(estimate_nhw_caide, estimate_asian_caide, 
                     estimate_black_caide, 
                     estimate_his_caide)
std_err_caide <- c(std_err_nhw_caide, std_err_asian_caide, std_err_black_caide, 
                    std_err_his_caide)

apoe <- c("NHW", "Asian", "Black", "Hispanic")
# apoe <- c("Low, e2+", "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+",
#          "High, e2+", "High, e3/e3", "High, e4+") 

data_caide <- tibble(apoe, estimate_caide, std_err_caide)
colnames(data_caide) <- c("apoe", "estimate_caide", "std_err_caide")

data_or_caide <- data_caide %>% mutate(
  lci_caide = estimate_caide - (std_err_caide * 1.96),
  hci_caide = estimate_caide + (std_err_caide * 1.96), 
  or_caide = exp(estimate_caide),
  or_lci_caide = exp(lci_caide),
  or_hci_caide = exp(hci_caide),
  colorx = case_when(
    (or_lci_caide > 1) | (or_hci_caide < 1) ~ "red"
  )
) 

# fix color error:
data_or_caide$colorx[is.na(data_or_caide$colorx)] <- "black"
data_or_e2$colorx[is.na(data_or_e2$colorx)] <- "black"
data_or_e4$colorx[is.na(data_or_e4$colorx)] <- "black"


var <- c("CAIDE", "CAIDE", "CAIDE", "CAIDE")
caide_facets1 = tibble(data_or_caide, var)

var <- c("e2+", "e2+", "e2+", "e2+")
e2_facets1 = tibble(data_or_e2, var)

var <- c("e4+", "e4+", "e4+", "e4+")
e4_facets1 = tibble(data_or_e4, var)


# rename and concatenate:
caide_facets <- dplyr::rename(caide_facets1, estimate = estimate_caide,
                               std_err = std_err_caide, lci = lci_caide,
                               hci = hci_caide, or = or_caide, race = apoe,
                               or_lci = or_lci_caide, or_hci = or_hci_caide)
e2_facets <- dplyr::rename(e2_facets1, estimate = estimate_e2,
                           std_err = std_err_e2, lci = lci_e2,
                           hci = hci_e2, or = or_e2, or_lci = or_lci_e2,
                           or_hci = or_hci_e2, race = apoe)
e4_facets <- dplyr::rename(e4_facets1, estimate = estimate_e4,
                           std_err = std_err_e4, lci = lci_e4,
                           hci = hci_e4, or = or_e4, or_lci = or_lci_e4,
                           or_hci = or_hci_e4, race = apoe)

data_facets <- rbind(e2_facets, e4_facets, caide_facets)

ggplot(data_facets, aes(x=or, y = factor(race, levels = c('NHW', 'Asian', 
                                                          'Black', 'Hispanic')), 
                        group=var, color = colorx)) +
  facet_grid(. ~ var, scales='free_x') +
  #labs(title = "CAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci, xmax = or_hci), 
                 linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c( "black", "red")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5)
  )

########################################################
# CAIDE: Slide 29
########################################################

caide_calc_nc <- read_csv("caide_calc_nc.csv") %>% 
  janitor::clean_names()

# apoe <- c("Low, e2+", "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+",
#          "High, e2+", "High, e3/e3", "High, e4+") 

# Re-do filters:
nhw_caide <-  dplyr::filter(caide_calc_nc, 
                            race == "Non-Hispanic White")
asian_caide <-  dplyr::filter(caide_calc_nc, 
                              race == "Asian")
black_caide <-  dplyr::filter(caide_calc_nc, 
                              race == "Black")
hispanic_caide <-  dplyr::filter(caide_calc_nc, 
                                 race == "Hispanic")

adni_caide_binary_nhw <- nhw_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_asian <- asian_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_black <- black_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_caide_binary_his <- hispanic_caide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ caide + apoe
adni_caide_binary_nhw$caide_apoe <- factor(adni_caide_binary_nhw$caide_apoe, 
                                     ordered = FALSE )
adni_caide_binary_nhw$caide_apoe = relevel(adni_caide_binary_nhw$caide_apoe, 
                                     ref = "Mid, e3/e3")
adni_caide_binary_asian$caide_apoe <- factor(adni_caide_binary_asian$caide_apoe, 
                                       ordered = FALSE )
adni_caide_binary_asian$caide_apoe = relevel(adni_caide_binary_asian$caide_apoe, 
                                       ref = "Mid, e3/e3")
adni_caide_binary_black$caide_apoe <- factor(adni_caide_binary_black$caide_apoe, 
                                       ordered = FALSE )
adni_caide_binary_black$caide_apoe = relevel(adni_caide_binary_black$caide_apoe, 
                                       ref = "Mid, e3/e3")
adni_caide_binary_his$caide_apoe <- factor(adni_caide_binary_his$caide_apoe, 
                                     ordered = FALSE )
adni_caide_binary_his$caide_apoe = relevel(adni_caide_binary_his$caide_apoe, 
                                     ref = "Mid, e3/e3")

# 1. logistic model: diagnosis ~ caide + apoe (without cholesterol)
model_caide_logistic_nhw <- glm(diag ~ caide_apoe, 
                                data=adni_caide_binary_nhw, 
                                family = binomial)
model_caide_logistic_asian <- glm(diag ~ caide_apoe, 
                                  data=adni_caide_binary_asian, 
                                  family = binomial)
model_caide_logistic_black <- glm(diag ~ caide_apoe, 
                                  data=adni_caide_binary_black, 
                                  family = binomial)
model_caide_logistic_his <- glm(diag ~ caide_apoe, 
                                data=adni_caide_binary_his, 
                                family = binomial)

# results give estimate and std.error for nhw:
results_nhw <- summary(model_caide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_high_e2 <- df_apoe_nhw$Estimate[2]
std_err_nhw_high_e2 <- df_apoe_nhw$`Std. Error`[2]
estimate_nhw_high_e3 <- df_apoe_nhw$Estimate[3]
std_err_nhw_high_e3 <- df_apoe_nhw$`Std. Error`[3]
estimate_nhw_high_e4 <- df_apoe_nhw$Estimate[4]
std_err_nhw_high_e4 <- df_apoe_nhw$`Std. Error`[4]
estimate_nhw_low_e2 <- df_apoe_nhw$Estimate[5]
std_err_nhw_low_e2 <- df_apoe_nhw$`Std. Error`[5]
estimate_nhw_low_e3 <- df_apoe_nhw$Estimate[6]
std_err_nhw_low_e3 <- df_apoe_nhw$`Std. Error`[6]
estimate_nhw_low_e4 <- df_apoe_nhw$Estimate[7]
std_err_nhw_low_e4 <- df_apoe_nhw$`Std. Error`[7]
estimate_nhw_mid_e2 <- df_apoe_nhw$Estimate[8]
std_err_nhw_mid_e2 <- df_apoe_nhw$`Std. Error`[8]
estimate_nhw_mid_e4 <- df_apoe_nhw$Estimate[9]
std_err_nhw_mid_e4 <- df_apoe_nhw$`Std. Error`[9]
estimate_nhw_mid_e3 <- 0.00
std_err_nhw_mid_e3 <- 0.00

# results give estimate and std.error for asian:
results_asian <- summary(model_caide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_high_e2 <- df_apoe_asian$Estimate[2]
std_err_asian_high_e2 <- df_apoe_asian$`Std. Error`[2]
estimate_asian_high_e3 <- df_apoe_asian$Estimate[3]
std_err_asian_high_e3 <- df_apoe_asian$`Std. Error`[3]
estimate_asian_high_e4 <- df_apoe_asian$Estimate[4]
std_err_asian_high_e4 <- df_apoe_asian$`Std. Error`[4]
estimate_asian_low_e2 <- df_apoe_asian$Estimate[5]
std_err_asian_low_e2 <- df_apoe_asian$`Std. Error`[5]
estimate_asian_low_e3 <- df_apoe_asian$Estimate[6]
std_err_asian_low_e3 <- df_apoe_asian$`Std. Error`[6]
estimate_asian_low_e4 <- df_apoe_asian$Estimate[7]
std_err_asian_low_e4 <- df_apoe_asian$`Std. Error`[7]
estimate_asian_mid_e2 <- df_apoe_asian$Estimate[8]
std_err_asian_mid_e2 <- df_apoe_asian$`Std. Error`[8]
estimate_asian_mid_e4 <- df_apoe_asian$Estimate[9]
std_err_asian_mid_e4 <- df_apoe_asian$`Std. Error`[9]
estimate_asian_mid_e3 <- 0.00
std_err_asian_mid_e3 <- 0.00

# results give estimate and std.error for black:
results_black <- summary(model_caide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_high_e2 <- df_apoe_black$Estimate[2]
std_err_black_high_e2 <- df_apoe_black$`Std. Error`[2]
estimate_black_high_e3 <- df_apoe_black$Estimate[3]
std_err_black_high_e3 <- df_apoe_black$`Std. Error`[3]
estimate_black_high_e4 <- df_apoe_black$Estimate[4]
std_err_black_high_e4 <- df_apoe_black$`Std. Error`[4]
estimate_black_low_e2 <- df_apoe_black$Estimate[5]
std_err_black_low_e2 <- df_apoe_black$`Std. Error`[5]
estimate_black_low_e3 <- df_apoe_black$Estimate[6]
std_err_black_low_e3 <- df_apoe_black$`Std. Error`[6]
estimate_black_low_e4 <- df_apoe_black$Estimate[7]
std_err_black_low_e4 <- df_apoe_black$`Std. Error`[7]
estimate_black_mid_e2 <- df_apoe_black$Estimate[8]
std_err_black_mid_e2 <- df_apoe_black$`Std. Error`[8]
estimate_black_mid_e4 <- df_apoe_black$Estimate[9]
std_err_black_mid_e4 <- df_apoe_black$`Std. Error`[9]
estimate_black_mid_e3 <- 0.00
std_err_black_mid_e3 <- 0.00

# results give estimate and std.error for hispanic:
results_his <- summary(model_caide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_his)
estimate_his_high_e2 <- df_apoe_his$Estimate[2]
std_err_his_high_e2 <- df_apoe_his$`Std. Error`[2]
estimate_his_high_e3 <- df_apoe_his$Estimate[3]
std_err_his_high_e3 <- df_apoe_his$`Std. Error`[3]
estimate_his_high_e4 <- df_apoe_his$Estimate[4]
std_err_his_high_e4 <- df_apoe_his$`Std. Error`[4]
estimate_his_low_e2 <- df_apoe_his$Estimate[5]
std_err_his_low_e2 <- df_apoe_his$`Std. Error`[5]
estimate_his_low_e3 <- df_apoe_his$Estimate[6]
std_err_his_low_e3 <- df_apoe_his$`Std. Error`[6]
estimate_his_low_e4 <- df_apoe_his$Estimate[7]
std_err_his_low_e4 <- df_apoe_his$`Std. Error`[7]
estimate_his_mid_e2 <- df_apoe_his$Estimate[8]
std_err_his_mid_e2 <- df_apoe_his$`Std. Error`[8]
estimate_his_mid_e4 <- df_apoe_his$Estimate[9]
std_err_his_mid_e4 <- df_apoe_his$`Std. Error`[9]
estimate_his_mid_e3 <- 0.00
std_err_his_mid_e3 <- 0.00

estimate_nhw <- c(estimate_nhw_low_e2, estimate_nhw_low_e3,
                    estimate_nhw_low_e4, estimate_nhw_mid_e2,
                    estimate_nhw_mid_e3, estimate_nhw_mid_e4,
                    estimate_nhw_high_e2, estimate_nhw_high_e3,
                    estimate_nhw_high_e4)

std_err_nhw <- c(std_err_nhw_low_e2, std_err_nhw_low_e3,
                   std_err_nhw_low_e4, std_err_nhw_mid_e2,
                   std_err_nhw_mid_e3, std_err_nhw_mid_e4,
                   std_err_nhw_high_e2, std_err_nhw_high_e3,
                   std_err_nhw_high_e4)

estimate_asian <- c(estimate_asian_low_e2, estimate_asian_low_e3,
                  estimate_asian_low_e4, estimate_asian_mid_e2,
                  estimate_asian_mid_e3, estimate_asian_mid_e4,
                  estimate_asian_high_e2, estimate_asian_high_e3,
                  estimate_asian_high_e4)

std_err_asian <- c(std_err_asian_low_e2, std_err_asian_low_e3,
                 std_err_asian_low_e4, std_err_asian_mid_e2,
                 std_err_asian_mid_e3, std_err_asian_mid_e4,
                 std_err_asian_high_e2, std_err_asian_high_e3,
                 std_err_asian_high_e4)

estimate_black <- c(estimate_black_low_e2, estimate_black_low_e3,
                    estimate_black_low_e4, estimate_black_mid_e2,
                    estimate_black_mid_e3, estimate_black_mid_e4,
                    estimate_black_high_e2, estimate_black_high_e3,
                    estimate_black_high_e4)

std_err_black <- c(std_err_black_low_e2, std_err_black_low_e3,
                   std_err_black_low_e4, std_err_black_mid_e2,
                   std_err_black_mid_e3, std_err_black_mid_e4,
                   std_err_black_high_e2, std_err_black_high_e3,
                   std_err_black_high_e4)

estimate_his <- c(estimate_his_low_e2, estimate_his_low_e3,
                    estimate_his_low_e4, estimate_his_mid_e2,
                    estimate_his_mid_e3, estimate_his_mid_e4,
                    estimate_his_high_e2, estimate_his_high_e3,
                    estimate_his_high_e4)

std_err_his <- c(std_err_his_low_e2, std_err_his_low_e3,
                   std_err_his_low_e4, std_err_his_mid_e2,
                   std_err_his_mid_e3, std_err_his_mid_e4,
                   std_err_his_high_e2, std_err_his_high_e3,
                   std_err_his_high_e4)

# apoe is the y-axis variables, e2 is the x-axis headers, so black, nhw, etc
apoe <- c("Low, e2+", "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e3/e3",
          "Mid, e4+", "High, e2+", "High, e3/e3", "High, e4+") 

data_nhw <- tibble(apoe, estimate_nhw, std_err_nhw)
data_asian <- tibble(apoe, estimate_asian, std_err_asian)
data_black <- tibble(apoe, estimate_black, std_err_black)
data_his <- tibble(apoe, estimate_his, std_err_his)

colnames(data_nhw) <- c("apoe", "estimate_nhw", "std_err_nhw")
colnames(data_asian) <- c("apoe", "estimate_asian", "std_err_asian")
colnames(data_black) <- c("apoe", "estimate_black", "std_err_black")
colnames(data_his) <- c("apoe", "estimate_his", "std_err_his")

data_or_nhw <- data_nhw %>% mutate(
  lci_nhw = estimate_nhw - (std_err_nhw * 1.96),
  hci_nhw = estimate_nhw + (std_err_nhw * 1.96), 
  or_nhw = exp(estimate_nhw),
  or_lci_nhw = exp(lci_nhw),
  or_hci_nhw = exp(hci_nhw),
  colorx = case_when(
    (or_lci_nhw > 1) | (or_hci_nhw < 1) ~ "red"
  )
)

data_or_asian <- data_asian %>% mutate(
  lci_asian = estimate_asian - (std_err_asian * 1.96),
  hci_asian = estimate_asian + (std_err_asian * 1.96), 
  or_asian = exp(estimate_asian),
  or_lci_asian = exp(lci_asian),
  or_hci_asian = exp(hci_asian),
  colorx = case_when(
    (or_lci_asian > 1) | (or_hci_asian < 1) ~ "red"
  )
)

data_or_black <- data_black %>% mutate(
  lci_black = estimate_black - (std_err_black * 1.96),
  hci_black = estimate_black + (std_err_black * 1.96), 
  or_black = exp(estimate_black),
  or_lci_black = exp(lci_black),
  or_hci_black = exp(hci_black),
  colorx = case_when(
    (or_lci_black > 1) | (or_hci_black < 1) ~ "red"
  )
)

data_or_his <- data_his %>% mutate(
  lci_his = estimate_his - (std_err_his * 1.96),
  hci_his = estimate_his + (std_err_his * 1.96), 
  or_his = exp(estimate_his),
  or_lci_his = exp(lci_his),
  or_hci_his = exp(hci_his),
  colorx = case_when(
    (or_lci_his > 1) | (or_hci_his < 1) ~ "red"
  )
)

# fix color error:
data_or_nhw$colorx[is.na(data_or_nhw$colorx)] <- "black"
data_or_asian$colorx[is.na(data_or_asian$colorx)] <- "black"
data_or_black$colorx[is.na(data_or_black$colorx)] <- "black"
data_or_his$colorx[is.na(data_or_his$colorx)] <- "black"

var <- c("NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW")
nhw_facets1 = tibble(data_or_nhw, var)

var <- c("Asian", "Asian", "Asian", "Asian", "Asian", "Asian", "Asian", 
         "Asian", "Asian")
asian_facets1 = tibble(data_or_asian, var)

var <- c("Black", "Black", "Black", "Black", "Black", "Black", "Black", 
         "Black", "Black")
black_facets1 = tibble(data_or_black, var)

var <- c("Hispanic", "Hispanic", "Hispanic", "Hispanic", "Hispanic", 
         "Hispanic", "Hispanic", "Hispanic", "Hispanic")
his_facets1 = tibble(data_or_his, var)

nhw_facets <- dplyr::rename(nhw_facets1, estimate = estimate_nhw,
                              std_err = std_err_nhw, lci = lci_nhw,
                              hci = hci_nhw, or = or_nhw, 
                              or_lci = or_lci_nhw, or_hci = or_hci_nhw)
asian_facets <- dplyr::rename(asian_facets1, estimate = estimate_asian,
                            std_err = std_err_asian, lci = lci_asian,
                            hci = hci_asian, or = or_asian, 
                            or_lci = or_lci_asian, or_hci = or_hci_asian)
black_facets <- dplyr::rename(black_facets1, estimate = estimate_black,
                              std_err = std_err_black, lci = lci_black,
                              hci = hci_black, or = or_black, 
                              or_lci = or_lci_black, or_hci = or_hci_black)
his_facets <- dplyr::rename(his_facets1, estimate = estimate_his,
                              std_err = std_err_his, lci = lci_his,
                              hci = hci_his, or = or_his, 
                              or_lci = or_lci_his, or_hci = or_hci_his)

data_facets <- rbind(nhw_facets, asian_facets, black_facets, his_facets)

ggplot(data_facets, aes(x=or, y = factor(apoe, levels = c("Low, e2+", 
                                                          "Low, e3/e3", 
                                                          "Low, e4+", 
                                                          "Mid, e2+", 
                                                          "Mid, e3/e3",
                                                          "Mid, e4+", 
                                                          "High, e2+", 
                                                          "High, e3/e3", 
                                                          "High, e4+")), 
                        group=var, color = colorx)) +
  facet_grid(. ~ var, scales='free_x') +
  #labs(title = "CAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci, xmax = or_hci), 
                 linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c( "black", "red")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5)
  )

########################################################
# mCAIDE: Slide 29
########################################################

mcaide_calc_nc <- read_csv("mcaide_calc_nc.csv") %>% 
  janitor::clean_names()

# apoe <- c("Low, e2+", "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e4+",
#          "High, e2+", "High, e3/e3", "High, e4+") 

# Re-do filters:
nhw_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                            race == "Non-Hispanic White")
asian_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                              race == "Asian")
black_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                              race == "Black")
hispanic_mcaide <-  dplyr::filter(mcaide_calc_nc, 
                                 race == "Hispanic")

adni_mcaide_binary_nhw <- nhw_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_asian <- asian_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_black <- black_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

adni_mcaide_binary_his <- hispanic_mcaide %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference for apoe:
# Model: diag ~ mcaide + apoe
adni_mcaide_binary_nhw$mcaide_apoe <- factor(adni_mcaide_binary_nhw$mcaide_apoe, 
                                           ordered = FALSE )
adni_mcaide_binary_nhw$mcaide_apoe = relevel(adni_mcaide_binary_nhw$mcaide_apoe, 
                                           ref = "Mid, e3/e3")
adni_mcaide_binary_asian$mcaide_apoe <- factor(adni_mcaide_binary_asian$mcaide_apoe, 
                                             ordered = FALSE )
adni_mcaide_binary_asian$mcaide_apoe = relevel(adni_mcaide_binary_asian$mcaide_apoe, 
                                             ref = "Mid, e3/e3")
adni_mcaide_binary_black$mcaide_apoe <- factor(adni_mcaide_binary_black$mcaide_apoe, 
                                             ordered = FALSE )
adni_mcaide_binary_black$mcaide_apoe = relevel(adni_mcaide_binary_black$mcaide_apoe, 
                                             ref = "Mid, e3/e3")
adni_mcaide_binary_his$mcaide_apoe <- factor(adni_mcaide_binary_his$mcaide_apoe, 
                                           ordered = FALSE )
adni_mcaide_binary_his$mcaide_apoe = relevel(adni_mcaide_binary_his$mcaide_apoe, 
                                           ref = "Mid, e3/e3")

# 1. logistic model: diagnosis ~ mcaide + apoe (without cholesterol)
model_mcaide_logistic_nhw <- glm(diag ~ mcaide_apoe, 
                                data=adni_mcaide_binary_nhw, 
                                family = binomial)
model_mcaide_logistic_asian <- glm(diag ~ mcaide_apoe, 
                                  data=adni_mcaide_binary_asian, 
                                  family = binomial)
model_mcaide_logistic_black <- glm(diag ~ mcaide_apoe, 
                                  data=adni_mcaide_binary_black, 
                                  family = binomial)
model_mcaide_logistic_his <- glm(diag ~ mcaide_apoe, 
                                data=adni_mcaide_binary_his, 
                                family = binomial)

# results give estimate and std.error for nhw:
results_nhw <- summary(model_mcaide_logistic_nhw)$coefficients
df_apoe_nhw <- as.data.frame(results_nhw)
estimate_nhw_high_e2 <- df_apoe_nhw$Estimate[2]
std_err_nhw_high_e2 <- df_apoe_nhw$`Std. Error`[2]
estimate_nhw_high_e3 <- df_apoe_nhw$Estimate[3]
std_err_nhw_high_e3 <- df_apoe_nhw$`Std. Error`[3]
estimate_nhw_high_e4 <- df_apoe_nhw$Estimate[4]
std_err_nhw_high_e4 <- df_apoe_nhw$`Std. Error`[4]
estimate_nhw_low_e2 <- df_apoe_nhw$Estimate[5]
std_err_nhw_low_e2 <- df_apoe_nhw$`Std. Error`[5]
estimate_nhw_low_e3 <- df_apoe_nhw$Estimate[6]
std_err_nhw_low_e3 <- df_apoe_nhw$`Std. Error`[6]
estimate_nhw_low_e4 <- df_apoe_nhw$Estimate[7]
std_err_nhw_low_e4 <- df_apoe_nhw$`Std. Error`[7]
estimate_nhw_mid_e2 <- df_apoe_nhw$Estimate[8]
std_err_nhw_mid_e2 <- df_apoe_nhw$`Std. Error`[8]
estimate_nhw_mid_e4 <- df_apoe_nhw$Estimate[9]
std_err_nhw_mid_e4 <- df_apoe_nhw$`Std. Error`[9]
estimate_nhw_mid_e3 <- 0.00
std_err_nhw_mid_e3 <- 0.00

# results give estimate and std.error for asian:
results_asian <- summary(model_mcaide_logistic_asian)$coefficients
df_apoe_asian <- as.data.frame(results_asian)
estimate_asian_high_e2 <- df_apoe_asian$Estimate[2]
std_err_asian_high_e2 <- df_apoe_asian$`Std. Error`[2]
estimate_asian_high_e3 <- df_apoe_asian$Estimate[3]
std_err_asian_high_e3 <- df_apoe_asian$`Std. Error`[3]
estimate_asian_high_e4 <- df_apoe_asian$Estimate[4]
std_err_asian_high_e4 <- df_apoe_asian$`Std. Error`[4]
estimate_asian_low_e2 <- df_apoe_asian$Estimate[5]
std_err_asian_low_e2 <- df_apoe_asian$`Std. Error`[5]
estimate_asian_low_e3 <- df_apoe_asian$Estimate[6]
std_err_asian_low_e3 <- df_apoe_asian$`Std. Error`[6]
estimate_asian_low_e4 <- df_apoe_asian$Estimate[7]
std_err_asian_low_e4 <- df_apoe_asian$`Std. Error`[7]
estimate_asian_mid_e2 <- df_apoe_asian$Estimate[8]
std_err_asian_mid_e2 <- df_apoe_asian$`Std. Error`[8]
estimate_asian_mid_e4 <- df_apoe_asian$Estimate[9]
std_err_asian_mid_e4 <- df_apoe_asian$`Std. Error`[9]
estimate_asian_mid_e3 <- 0.00
std_err_asian_mid_e3 <- 0.00

# results give estimate and std.error for black:
results_black <- summary(model_mcaide_logistic_black)$coefficients
df_apoe_black <- as.data.frame(results_black)
estimate_black_high_e2 <- df_apoe_black$Estimate[2]
std_err_black_high_e2 <- df_apoe_black$`Std. Error`[2]
estimate_black_high_e3 <- df_apoe_black$Estimate[3]
std_err_black_high_e3 <- df_apoe_black$`Std. Error`[3]
estimate_black_high_e4 <- df_apoe_black$Estimate[4]
std_err_black_high_e4 <- df_apoe_black$`Std. Error`[4]
estimate_black_low_e2 <- df_apoe_black$Estimate[5]
std_err_black_low_e2 <- df_apoe_black$`Std. Error`[5]
estimate_black_low_e3 <- df_apoe_black$Estimate[6]
std_err_black_low_e3 <- df_apoe_black$`Std. Error`[6]
estimate_black_low_e4 <- df_apoe_black$Estimate[7]
std_err_black_low_e4 <- df_apoe_black$`Std. Error`[7]
estimate_black_mid_e2 <- df_apoe_black$Estimate[8]
std_err_black_mid_e2 <- df_apoe_black$`Std. Error`[8]
estimate_black_mid_e4 <- df_apoe_black$Estimate[9]
std_err_black_mid_e4 <- df_apoe_black$`Std. Error`[9]
estimate_black_mid_e3 <- 0.00
std_err_black_mid_e3 <- 0.00

# results give estimate and std.error for hispanic:
results_his <- summary(model_mcaide_logistic_his)$coefficients
df_apoe_his <- as.data.frame(results_his)
estimate_his_high_e2 <- df_apoe_his$Estimate[2]
std_err_his_high_e2 <- df_apoe_his$`Std. Error`[2]
estimate_his_high_e3 <- df_apoe_his$Estimate[3]
std_err_his_high_e3 <- df_apoe_his$`Std. Error`[3]
estimate_his_high_e4 <- df_apoe_his$Estimate[4]
std_err_his_high_e4 <- df_apoe_his$`Std. Error`[4]
estimate_his_low_e2 <- df_apoe_his$Estimate[5]
std_err_his_low_e2 <- df_apoe_his$`Std. Error`[5]
estimate_his_low_e3 <- df_apoe_his$Estimate[6]
std_err_his_low_e3 <- df_apoe_his$`Std. Error`[6]
estimate_his_low_e4 <- df_apoe_his$Estimate[7]
std_err_his_low_e4 <- df_apoe_his$`Std. Error`[7]
estimate_his_mid_e2 <- df_apoe_his$Estimate[8]
std_err_his_mid_e2 <- df_apoe_his$`Std. Error`[8]
estimate_his_mid_e4 <- df_apoe_his$Estimate[9]
std_err_his_mid_e4 <- df_apoe_his$`Std. Error`[9]
estimate_his_mid_e3 <- 0.00
std_err_his_mid_e3 <- 0.00

estimate_nhw <- c(estimate_nhw_low_e2, estimate_nhw_low_e3,
                  estimate_nhw_low_e4, estimate_nhw_mid_e2,
                  estimate_nhw_mid_e3, estimate_nhw_mid_e4,
                  estimate_nhw_high_e2, estimate_nhw_high_e3,
                  estimate_nhw_high_e4)

std_err_nhw <- c(std_err_nhw_low_e2, std_err_nhw_low_e3,
                 std_err_nhw_low_e4, std_err_nhw_mid_e2,
                 std_err_nhw_mid_e3, std_err_nhw_mid_e4,
                 std_err_nhw_high_e2, std_err_nhw_high_e3,
                 std_err_nhw_high_e4)

estimate_asian <- c(estimate_asian_low_e2, estimate_asian_low_e3,
                    estimate_asian_low_e4, estimate_asian_mid_e2,
                    estimate_asian_mid_e3, estimate_asian_mid_e4,
                    estimate_asian_high_e2, estimate_asian_high_e3,
                    estimate_asian_high_e4)

std_err_asian <- c(std_err_asian_low_e2, std_err_asian_low_e3,
                   std_err_asian_low_e4, std_err_asian_mid_e2,
                   std_err_asian_mid_e3, std_err_asian_mid_e4,
                   std_err_asian_high_e2, std_err_asian_high_e3,
                   std_err_asian_high_e4)

estimate_black <- c(estimate_black_low_e2, estimate_black_low_e3,
                    estimate_black_low_e4, estimate_black_mid_e2,
                    estimate_black_mid_e3, estimate_black_mid_e4,
                    estimate_black_high_e2, estimate_black_high_e3,
                    estimate_black_high_e4)

std_err_black <- c(std_err_black_low_e2, std_err_black_low_e3,
                   std_err_black_low_e4, std_err_black_mid_e2,
                   std_err_black_mid_e3, std_err_black_mid_e4,
                   std_err_black_high_e2, std_err_black_high_e3,
                   std_err_black_high_e4)

estimate_his <- c(estimate_his_low_e2, estimate_his_low_e3,
                  estimate_his_low_e4, estimate_his_mid_e2,
                  estimate_his_mid_e3, estimate_his_mid_e4,
                  estimate_his_high_e2, estimate_his_high_e3,
                  estimate_his_high_e4)

std_err_his <- c(std_err_his_low_e2, std_err_his_low_e3,
                 std_err_his_low_e4, std_err_his_mid_e2,
                 std_err_his_mid_e3, std_err_his_mid_e4,
                 std_err_his_high_e2, std_err_his_high_e3,
                 std_err_his_high_e4)

# apoe is the y-axis variables, e2 is the x-axis headers, so black, nhw, etc
apoe <- c("Low, e2+", "Low, e3/e3", "Low, e4+", "Mid, e2+", "Mid, e3/e3",
          "Mid, e4+", "High, e2+", "High, e3/e3", "High, e4+") 

data_nhw <- tibble(apoe, estimate_nhw, std_err_nhw)
data_asian <- tibble(apoe, estimate_asian, std_err_asian)
data_black <- tibble(apoe, estimate_black, std_err_black)
data_his <- tibble(apoe, estimate_his, std_err_his)

colnames(data_nhw) <- c("apoe", "estimate_nhw", "std_err_nhw")
colnames(data_asian) <- c("apoe", "estimate_asian", "std_err_asian")
colnames(data_black) <- c("apoe", "estimate_black", "std_err_black")
colnames(data_his) <- c("apoe", "estimate_his", "std_err_his")

data_or_nhw <- data_nhw %>% mutate(
  lci_nhw = estimate_nhw - (std_err_nhw * 1.96),
  hci_nhw = estimate_nhw + (std_err_nhw * 1.96), 
  or_nhw = exp(estimate_nhw),
  or_lci_nhw = exp(lci_nhw),
  or_hci_nhw = exp(hci_nhw),
  colorx = case_when(
    (or_lci_nhw > 1) | (or_hci_nhw < 1) ~ "red"
  )
)

data_or_asian <- data_asian %>% mutate(
  lci_asian = estimate_asian - (std_err_asian * 1.96),
  hci_asian = estimate_asian + (std_err_asian * 1.96), 
  or_asian = exp(estimate_asian),
  or_lci_asian = exp(lci_asian),
  or_hci_asian = exp(hci_asian),
  colorx = case_when(
    (or_lci_asian > 1) | (or_hci_asian < 1) ~ "red"
  )
)

data_or_black <- data_black %>% mutate(
  lci_black = estimate_black - (std_err_black * 1.96),
  hci_black = estimate_black + (std_err_black * 1.96), 
  or_black = exp(estimate_black),
  or_lci_black = exp(lci_black),
  or_hci_black = exp(hci_black),
  colorx = case_when(
    (or_lci_black > 1) | (or_hci_black < 1) ~ "red"
  )
)

data_or_his <- data_his %>% mutate(
  lci_his = estimate_his - (std_err_his * 1.96),
  hci_his = estimate_his + (std_err_his * 1.96), 
  or_his = exp(estimate_his),
  or_lci_his = exp(lci_his),
  or_hci_his = exp(hci_his),
  colorx = case_when(
    (or_lci_his > 1) | (or_hci_his < 1) ~ "red"
  )
)

# fix color error:
data_or_nhw$colorx[is.na(data_or_nhw$colorx)] <- "black"
data_or_asian$colorx[is.na(data_or_asian$colorx)] <- "black"
data_or_black$colorx[is.na(data_or_black$colorx)] <- "black"
data_or_his$colorx[is.na(data_or_his$colorx)] <- "black"

var <- c("NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW", "NHW")
nhw_facets1 = tibble(data_or_nhw, var)

var <- c("Asian", "Asian", "Asian", "Asian", "Asian", "Asian", "Asian", 
         "Asian", "Asian")
asian_facets1 = tibble(data_or_asian, var)

var <- c("Black", "Black", "Black", "Black", "Black", "Black", "Black", 
         "Black", "Black")
black_facets1 = tibble(data_or_black, var)

var <- c("Hispanic", "Hispanic", "Hispanic", "Hispanic", "Hispanic", 
         "Hispanic", "Hispanic", "Hispanic", "Hispanic")
his_facets1 = tibble(data_or_his, var)

nhw_facets <- dplyr::rename(nhw_facets1, estimate = estimate_nhw,
                            std_err = std_err_nhw, lci = lci_nhw,
                            hci = hci_nhw, or = or_nhw, 
                            or_lci = or_lci_nhw, or_hci = or_hci_nhw)
asian_facets <- dplyr::rename(asian_facets1, estimate = estimate_asian,
                              std_err = std_err_asian, lci = lci_asian,
                              hci = hci_asian, or = or_asian, 
                              or_lci = or_lci_asian, or_hci = or_hci_asian)
black_facets <- dplyr::rename(black_facets1, estimate = estimate_black,
                              std_err = std_err_black, lci = lci_black,
                              hci = hci_black, or = or_black, 
                              or_lci = or_lci_black, or_hci = or_hci_black)
his_facets <- dplyr::rename(his_facets1, estimate = estimate_his,
                            std_err = std_err_his, lci = lci_his,
                            hci = hci_his, or = or_his, 
                            or_lci = or_lci_his, or_hci = or_hci_his)

data_facets <- rbind(nhw_facets, asian_facets, black_facets, his_facets)

ggplot(data_facets, aes(x=or, y = factor(apoe, levels = c("Low, e2+", 
                                                          "Low, e3/e3", 
                                                          "Low, e4+", 
                                                          "Mid, e2+", 
                                                          "Mid, e3/e3",
                                                          "Mid, e4+", 
                                                          "High, e2+", 
                                                          "High, e3/e3", 
                                                          "High, e4+")), 
                        group=var, color = colorx)) +
  facet_grid(. ~ var, scales='free_x') +
  #labs(title = "CAIDE") +
  geom_vline(xintercept = 1, linetype = 2) + 
  geom_point() +
  geom_linerange(aes(xmin = or_lci, xmax = or_hci), 
                 linewidth = 0.75) +
  scale_y_discrete(" ") +
  scale_x_continuous("Odds Ratio", labels = 
                       scales::number_format(accuracy = 0.1)) + 
  # scale_x_continuous("Odds Ratio", labels = 
  #                      scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c( "black", "red")) +
  theme_bw() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0),
    # axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=0),
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(face = 'bold'),
    plot.title = element_text(hjust = 0.5)
  )







