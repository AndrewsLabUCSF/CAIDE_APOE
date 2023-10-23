# setwd("/wynton/group/andrews/users/evelasquez/CAIDE_APOE/scripts/ADNI")
setwd("/Users/elinorvelasquez/Desktop/ADNI/")
library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("readxl")
library("readr")
# library("webshot2")
library("ggplot2")
library("gtsummary")
library("broom")
library("MASS")

####################################################
# Construct APOE
####################################################
apoe_raw <- read_csv("APOERES_28Sep2023.csv") %>% 
  janitor::clean_names()

# APOE redefined:
apoe_apgen <- apoe_raw %>%
  unite(apoe_geno, c(apgen1, apgen2), sep = "") %>%
  mutate(
    apoe = case_when(
      apoe_geno %in% c("22", "23", "32") ~ "e2+",
      apoe_geno %in% c("33") ~ "e3/e3",
      apoe_geno %in% c("34", "24", "44", "43", "42") ~ "e4+",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(c(rid, apoe))

apoe_apgen %>% count(apoe) # NAs 0
apoe_apgen %>% count(apoe_geno) #NAs 0
apoe_raw %>% count(apgen1) #NAs 0
apoe_raw %>% count(apgen2) #NAs 0
janitor::get_dupes(apoe_apgen) # No dupes

#####################################################
# compute caide_basic (adds "race")
#####################################################

merge_raw <- read_csv("~/Desktop/ADNI/ADNIMERGE_09Aug2023.csv") %>% 
  janitor::clean_names()

#### Add race:
caide_basic <- merge_raw %>%
  mutate(
    race = case_when(
      # ptraccat != "Asian" & ptethcat != "Hisp/Latino" & ptraccat != "Black" 
      #  & ptraccat != "White" ~ "Other",
      ptraccat == "More than one" ~ "Other",
      ptraccat == "Unknown" ~ "Other",
      ptraccat == "Am Indian/Alaskan" ~ "Other",
      ptraccat == "Hawaiian/Other PI" ~ "Other",
      ptraccat == "White" & ptethcat != "Hisp/Latino" ~ "NHW",
      ptraccat == "White" & ptethcat == "Hisp/Latino" ~ "Hispanic",
      ptethcat == "Unknown" ~ NA_character_,
      ptraccat == "Black" ~ "Black",
      ptraccat == "Asian" ~ "Asian",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(viscode == "bl") %>%
  dplyr::filter(dx_bl != "NA") %>% # Unknowns 11 sit in dx_bl = "NA"
  dplyr::filter(race != "Other") %>%
  dplyr::select(rid, viscode, age, ptgender, cdrsb,
                                  pteducat, race, dx_bl) %>%
  mutate(
    viscode = as.factor(viscode),
    sex_factor = as.factor(ptgender)
  ) %>%
  dplyr::select(-c(ptgender)) %>%
  dplyr::filter(viscode == 'bl') %>%
  dplyr::select(-c(viscode))

#janitor::get_dupes(caide_basic) # No dupes

###################################################################
# COMPUTE CHOLESTEROL
###################################################################

chol_rct20_raw <- read_csv("~/Desktop/ADNI/LABDATA_27Sep2023.csv") %>% 
  janitor::clean_names()

chol_rct20 <- chol_rct20_raw %>% dplyr::filter(viscode2 == "sc" | 
                        viscode2 == "bl") %>%
                        mutate(
                          rct20 = ifelse(rct20 == -1 | 
                                         rct20 == "ADC05" |
                                         rct20 == "ADC259" |
                                         rct20 == "SCC09" |
                                         rct20 == "SCC10" |
                                         rct20 == "SCC142" |
                                         rct20 == "SCC18",
                                         NA_character_, rct20)
                        ) %>%
                        dplyr::select(rid, examdate, rct20) %>%
                        group_by(rid) %>% arrange(rid, examdate) %>%
                        slice_min(examdate) %>% ungroup()

##str(chol_rct20) # 1920 rows

##chol_rct20a %>% dplyr::filter(rid==4848)
# rid  rct20 viscode2 examdate  
# 4848 156   sc       2012-07-05
# 4848 156   sc       2012-10-17

# Filter duplicates
##janitor::get_dupes(chol_rct20_raw) # no dupes
##janitor::get_dupes(chol_rct20a) # dupes
# rid = 4848 is duplicated --fixed: see above code

#####################################################
# COMPUTE BMI & SYSBP
#####################################################
vitals_raw <- read_csv("VITALS_14Aug2023.csv") %>% 
  janitor::clean_names() 

pre_bmi_weight <- vitals_raw %>% 
                      dplyr::filter(viscode2 == "bl") %>%
                      dplyr::select(rid, vsweight, vswtunit, vsbpsys)

pre_bmi_height <- vitals_raw %>% 
                      dplyr::filter(viscode2 == "sc") %>%
                      dplyr::select(rid, vsheight, vshtunit) 

pre_bmi_height %>% count(is.na(vsheight))
                                                     
pre_bmi <- pre_bmi_weight %>% left_join(pre_bmi_height, 
                                        by=c('rid'))

bmi_computed <- pre_bmi %>%
  mutate(
    weight_in_kg = case_when(vswtunit == 1 ~ vsweight*0.453592, # 1 is lbs
                             vswtunit == 2 ~ vsweight,
                             vswtunit == -1 | vswtunit == -4 ~ NA_real_,
                             vsweight == -1 | vsweight == -4 ~ NA_real_,
                             TRUE ~ NA_real_
    ),
    height_in_cm = case_when(vshtunit == 1 ~ vsheight*2.54, # 1 is inches
                             vshtunit == 2 ~ vsheight,
                             vshtunit == -1 | vshtunit == -4 ~ NA_real_,
                             vsheight == -1 | vsheight == -4 ~ NA_real_,
                             TRUE ~ NA_real_
    ),
    bmi = case_when(
      height_in_cm != "NA" & weight_in_kg != "NA" ~ 
        weight_in_kg/((height_in_cm/100)**2), # From CDC
      height_in_cm == "NA" | weight_in_kg == "NA" ~ NA_real_
    )
  ) %>%
  dplyr::select(rid, bmi, vsbpsys)

janitor::get_dupes(bmi_computed) # No dupes
bmi_computed %>% count(is.na(bmi)) # NAs 36 
########################################################
# Left joins: left_joined_data
########################################################

left_joined_data <- caide_basic %>% 
                        left_join(chol_rct20, by=c('rid')) %>%
                        left_join(bmi_computed, by=c('rid')) %>%
                        left_join(apoe_apgen, by=c('rid')) %>%
                        group_by(rid) %>% 
                        arrange(rid, examdate) %>%
                        slice_min(examdate) %>% ungroup()

janitor::get_dupes(left_joined_data) # no dupes
str(left_joined_data) # 2369 rows
left_joined_data %>% count(is.na(apoe)) # NAs 195

###################################################
# Imputing variables
# Case 3: Chol is imputed
# Case 1: Chol is not imputed; caide doesn't include chol: 
# Case 2: Chol is not imputed; caide includes chol (with NAs)
###################################################

# prep for imputation: Case 3
impute <- left_joined_data %>%
  ungroup() %>%
  dplyr::select(-c(rid, race, dx_bl, examdate)) %>%
  mutate(
    apoe_factor = as.factor(apoe),
    rct20_num = as.numeric(rct20)
  ) %>%
  dplyr::select(-c(apoe, rct20))

str(impute) #imputing cholesterol (rct20_num)
# age cdrsb pteducat sex_factor bmi vsbpsys apoe_factor rct20_num

not_imputed <- dplyr::select(left_joined_data, 
                        'rid', 'race', 'dx_bl', 'examdate', 'rct20')

janitor::get_dupes(left_joined_data) # no dupes
janitor::get_dupes(not_imputed) # no dupes
janitor::get_dupes(impute) # no dupes

# imputing...CASE 3: IMPUTED CHOL
impute_df <- as.data.frame(impute)
no_missing_values <- missForest(impute_df)
impute_nmv <- no_missing_values$ximp

# rejoin impute and not_imputed: Case 3
complete_vars <- cbind(not_imputed, impute_nmv) %>%
                        tibble() %>%
                        mutate(dx_bl = as.factor(dx_bl))

# replace categorical vars in imputed dataset: apoe, sex
# capture original apoe & sex
orig_apoe_sex <- dplyr::select(left_joined_data, 'rid', 'apoe', 
                                    'sex_factor', 'examdate') %>%
                      mutate(
                        apoe = as.factor(apoe),
                        sex_factor = as.factor(sex_factor)
                      ) 

# Filter duplicates
janitor::get_dupes(left_joined_data) # No dups
janitor::get_dupes(orig_apoe_sex) # no dupes

# subtract apoe & sex_factor from imputed dataset
minor <- complete_vars %>% 
  dplyr::select(-c(apoe_factor, sex_factor))

# add the unimputed apoe & sex to "minor"
full_vars <- minor %>%
  left_join(orig_apoe_sex, by=c('rid', 'examdate')) %>%
  group_by(rid) %>% 
  arrange(rid, examdate) %>%
  slice_min(examdate) %>%
  ungroup() %>%
  mutate(
    rct20_no_imp = as.numeric(rct20)
  ) %>%
  dplyr::select(-c('rct20'))
  

str(full_vars) # rct20_num (imputed chol); rct20_no_imp (not imputed)

###############################################################
# Case 3: Caide calc w/ cholesterol
###############################################################

caide_calc <- full_vars %>% 
  dplyr::rename(i_age = age, i_education =
                  pteducat, i_bmi = bmi, 
                i_systbp = vsbpsys, gender = sex_factor) %>%
  mutate(
    hypchol3 = case_when(     
        rct20_num >= 240 ~ 1,
        rct20_num < 240 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
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
      gender == "Male" ~ 1,
      gender == "Female" ~ 0,
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
    caide_chol3 = case_when(
      hypchol3 == 0 ~ 0, 
      hypchol3 == 1 ~ 2, 
      TRUE ~ NA_real_
    ),
  ) %>%
  rowwise() %>%
  mutate(
    caide3 = sum(caide_age, caide_educ, caide_sex, caide_obesity, 
                 caide_sbp, caide_chol3, na.rm = F),
    caide_missing3 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                         is.na(i_bmi), is.na(i_systbp))
  ) %>%
  ungroup() %>%
  mutate(
    z_caide3 = scale(caide3)[,1],
    caide_cat3 = case_when(
      between(caide3, 0, (mean(caide3, na.rm = T) - sd(caide3, na.rm = T))) 
      ~ 'Low',
      between(caide3, (mean(caide3, na.rm = T) - sd(caide3, na.rm = T)), 
              (mean(caide3, na.rm = T) + sd(caide3, na.rm = T))) ~ 'Mid',
      between(caide3, (mean(caide3, na.rm = T) + sd(caide3, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    caide_cat3 = fct_relevel(caide_cat3, 'Low', 'Mid',  'High'),
    caide_apoe3 = glue("{caide_cat3}, {apoe}"),
    caideXapoe3 = case_when(
      caide_apoe3 == "Low, NA" | caide_apoe3 == "Mid, NA" | 
          caide_apoe3 == "High, NA" ~ NA_character_,
      caide_apoe3 == "Low, e2+" ~ "Low, e2+",
      caide_apoe3 == "Low, e3/e3" ~ "Low, e3/e3",
      caide_apoe3 == "Low, e4+" ~ "Low, e4+",
      caide_apoe3 == "Mid, e2+" ~ "Mid, e2+",
      caide_apoe3 == "Mid, e3/e3" ~ "Mid, e3/e3",
      caide_apoe3 == "Mid, e4+" ~ "Mid, e4+",
      caide_apoe3 == "High, e2+" ~ "High, e2+",
      caide_apoe3 == "High, e3/e3" ~ "High, e3/e3",
      caide_apoe3 == "High, e4+" ~ "High, e4+"
    ),
    caide_apoe3 = fct_relevel(caideXapoe3, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                              "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                              "High, e4+")
  ) %>%
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
     gender == "Male" ~ 1,
     gender == "Female" ~ 0,
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
   ),
   mcaide_chol3 = case_when(
     hypchol3 == 0 ~ 0, 
     hypchol3 == 1 ~ 2, 
     TRUE ~ NA_real_
   )
 ) %>%
rowwise() %>%
  mutate(
    mcaide3 = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, 
                  mcaide_sbp, mcaide_chol3, na.rm = F),
    mcaide_missing3 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                          is.na(i_bmi), is.na(i_systbp), is.na(rct20_num))
  ) %>%                                                 
  ungroup() %>%
  mutate(
    z_mcaide3 = scale(mcaide3)[,1],
    mcaide_cat3 = case_when(
      between(mcaide3, 0, (mean(mcaide3, na.rm = T) - sd(mcaide3, na.rm = T))) 
      ~ 'Low',
      between(mcaide3, (mean(mcaide3, na.rm = T) - sd(mcaide3, na.rm = T)), 
              (mean(mcaide3, na.rm = T) + sd(mcaide3, na.rm = T))) ~ 'Mid',
      between(mcaide3, (mean(mcaide3, na.rm = T) + sd(mcaide3, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    mcaide_cat3 = fct_relevel(mcaide_cat3, 'Low', "Mid",  'High'),
    mcaide_apoe3 = glue("{mcaide_cat3}, {apoe}"),
    mcaideXapoe3 = case_when(
      mcaide_apoe3 == "Low, NA" | mcaide_apoe3 == "Mid, NA" | 
        mcaide_apoe3 == "High, NA" ~ NA_character_,
      mcaide_apoe3 == "Low, e2+" ~ "Low, e2+",
      mcaide_apoe3 == "Low, e3/e3" ~ "Low, e3/e3",
      mcaide_apoe3 == "Low, e4+" ~ "Low, e4+",
      mcaide_apoe3 == "Mid, e2+" ~ "Mid, e2+",
      mcaide_apoe3 == "Mid, e3/e3" ~ "Mid, e3/e3",
      mcaide_apoe3 == "Mid, e4+" ~ "Mid, e4+",
      mcaide_apoe3 == "High, e2+" ~ "High, e2+",
      mcaide_apoe3 == "High, e3/e3" ~ "High, e3/e3",
      mcaide_apoe3 == "High, e4+" ~ "High, e4+"
    ),
    mcaide_apoe3 = fct_relevel(mcaideXapoe3, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                               "Low, e4+",
                               "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                               "High, e4+")
  ) %>%
mutate(
  hypchol1 = case_when(
    rct20_no_imp >= 240 ~ 1,
    rct20_no_imp < 240 ~ 0,
    TRUE ~ NA_real_
  ),
  caide_chol1 = case_when(
    hypchol1 == 0 ~ 0, 
    hypchol1 == 1 ~ 2, 
    TRUE ~ NA_real_
  ),
) %>%
  rowwise() %>%
  mutate(
    caide1 = sum(caide_age, caide_educ, caide_sex, caide_obesity, 
                 caide_sbp, na.rm = F),
    caide_missing1 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                        is.na(i_bmi), is.na(i_systbp))
  ) %>%
  ungroup() %>%
  mutate(
    z_caide1 = scale(caide1)[,1],
    caide_cat1 = case_when(
      between(caide1, 0, (mean(caide1, na.rm = T) - sd(caide1, na.rm = T))) 
      ~ 'Low',
      between(caide1, (mean(caide1, na.rm = T) - sd(caide1, na.rm = T)), 
              (mean(caide1, na.rm = T) + sd(caide1, na.rm = T))) ~ 'Mid',
      between(caide1, (mean(caide1, na.rm = T) + sd(caide1, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    caide_cat1 = fct_relevel(caide_cat1, 'Low', 'Mid',  'High'),
    caide_apoe1 = glue("{caide_cat1}, {apoe}"),
    caideXapoe1 = case_when(
      caide_apoe1 == "Low, NA" | caide_apoe1 == "Mid, NA" | 
        caide_apoe1 == "High, NA" ~ NA_character_,
      caide_apoe1 == "Low, e2+" ~ "Low, e2+",
      caide_apoe1 == "Low, e3/e3" ~ "Low, e3/e3",
      caide_apoe1 == "Low, e4+" ~ "Low, e4+",
      caide_apoe1 == "Mid, e2+" ~ "Mid, e2+",
      caide_apoe1 == "Mid, e3/e3" ~ "Mid, e3/e3",
      caide_apoe1 == "Mid, e4+" ~ "Mid, e4+",
      caide_apoe1 == "High, e2+" ~ "High, e2+",
      caide_apoe1 == "High, e3/e3" ~ "High, e3/e3",
      caide_apoe1 == "High, e4+" ~ "High, e4+"
    ),
    caide_apoe1 = fct_relevel(caideXapoe1, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                             "Low, e4+",
                             "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                             "High, e4+")
  ) %>%
  rowwise() %>%
  mutate(
    mcaide1 = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, 
                 mcaide_sbp, na.rm = F),
    mcaide_missing1 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                         is.na(i_bmi), is.na(i_systbp))
  ) %>%
  ungroup() %>%
  mutate(
    z_mcaide1 = scale(mcaide1)[,1],
    mcaide_cat1 = case_when(
      between(mcaide1, 0, (mean(mcaide1, na.rm = T) - sd(mcaide1, na.rm = T))) 
      ~ 'Low',
      between(mcaide1, (mean(mcaide1, na.rm = T) - sd(mcaide1, na.rm = T)), 
              (mean(mcaide1, na.rm = T) + sd(mcaide1, na.rm = T))) ~ 'Mid',
      between(mcaide1, (mean(mcaide1, na.rm = T) + sd(mcaide1, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    mcaide_cat1 = fct_relevel(mcaide_cat1, 'Low', "Mid",  'High'),
    mcaide_apoe1 = glue("{mcaide_cat1}, {apoe}"),
    mcaideXapoe1 = case_when(
      mcaide_apoe1 == "Low, NA" | mcaide_apoe1 == "Mid, NA" | 
        mcaide_apoe1 == "High, NA" ~ NA_character_,
      mcaide_apoe1 == "Low, e2+" ~ "Low, e2+",
      mcaide_apoe1 == "Low, e3/e3" ~ "Low, e3/e3",
      mcaide_apoe1 == "Low, e4+" ~ "Low, e4+",
      mcaide_apoe1 == "Mid, e2+" ~ "Mid, e2+",
      mcaide_apoe1 == "Mid, e3/e3" ~ "Mid, e3/e3",
      mcaide_apoe1 == "Mid, e4+" ~ "Mid, e4+",
      mcaide_apoe1 == "High, e2+" ~ "High, e2+",
      mcaide_apoe1 == "High, e3/e3" ~ "High, e3/e3",
      mcaide_apoe1 == "High, e4+" ~ "High, e4+"
    ),
    mcaide_apoe1 = fct_relevel(mcaideXapoe1, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                              "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                              "High, e4+")
  ) %>%
#########################################################
# Case 2: Caide calc w/ raw chol (includes missing values)
#########################################################
  mutate(
    hypchol2 = case_when(
      rct20_no_imp >= 240 ~ 1,
      rct20_no_imp < 240 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    caide_chol2 = case_when(
      hypchol2 == 0 ~ 0, 
      hypchol2 == 1 ~ 2, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    caide2 = sum(caide_age, caide_educ, caide_sex, caide_obesity, 
                 caide_sbp, caide_chol2, na.rm = F),
    caide_missing2 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                        is.na(i_bmi), is.na(i_systbp), is.na(rct20_no_imp))
  ) %>%
  ungroup() %>%
  mutate(
    z_caide2 = scale(caide2)[,1],
    caide_cat2 = case_when(
      between(caide2, 0, (mean(caide2, na.rm = T) - sd(caide2, na.rm = T))) 
      ~ 'Low',
      between(caide2, (mean(caide2, na.rm = T) - sd(caide2, na.rm = T)), 
              (mean(caide2, na.rm = T) + sd(caide2, na.rm = T))) ~ 'Mid',
      between(caide2, (mean(caide2, na.rm = T) + sd(caide2, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    caide_cat2 = fct_relevel(caide_cat2, 'Low', 'Mid',  'High'),
    caide_apoe2 = glue("{caide_cat2}, {apoe}"),
    caideXapoe2 = case_when(
      caide_apoe2 == "Low, NA" | caide_apoe2 == "Mid, NA" | 
        caide_apoe2 == "High, NA" ~ NA_character_,
      caide_apoe2 == "Low, e2+" ~ "Low, e2+",
      caide_apoe2 == "Low, e3/e3" ~ "Low, e3/e3",
      caide_apoe2 == "Low, e4+" ~ "Low, e4+",
      caide_apoe2 == "Mid, e2+" ~ "Mid, e2+",
      caide_apoe2 == "Mid, e3/e3" ~ "Mid, e3/e3",
      caide_apoe2 == "Mid, e4+" ~ "Mid, e4+",
      caide_apoe2 == "High, e2+" ~ "High, e2+",
      caide_apoe2 == "High, e3/e3" ~ "High, e3/e3",
      caide_apoe2 == "High, e4+" ~ "High, e4+"
    ),
    caide_apoe2 = fct_relevel(caideXapoe2, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                             "Low, e4+",
                             "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                             "High, e4+")
  ) %>%
  mutate(
    mcaide_chol2 = case_when(
      hypchol2 == 0 ~ 0, 
      hypchol2 == 1 ~ 2, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    mcaide2 = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, 
                 mcaide_sbp, mcaide_chol2, na.rm = F),
    mcaide_missing2 = sum(is.na(i_age), is.na(i_education), is.na(gender), 
                         is.na(i_bmi), is.na(i_systbp), is.na(rct20_no_imp))
  ) %>%
  ungroup() %>%
  mutate(
    z_mcaide2 = scale(mcaide2)[,1],
    mcaide_cat2 = case_when(
      between(mcaide2, 0, (mean(mcaide2, na.rm = T) - sd(mcaide2, na.rm = T))) 
      ~ 'Low',
      between(mcaide2, (mean(mcaide2, na.rm = T) - sd(mcaide2, na.rm = T)), 
              (mean(mcaide2, na.rm = T) + sd(mcaide2, na.rm = T))) ~ 'Mid',
      between(mcaide2, (mean(mcaide2, na.rm = T) + sd(mcaide2, na.rm = T)), 14) ~ 
        'High',
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    mcaide_cat2 = fct_relevel(mcaide_cat2, 'Low', "Mid",  'High'),
    mcaide_apoe2 = glue("{mcaide_cat2}, {apoe}"),
    mcaideXapoe2 = case_when(
      mcaide_apoe2 == "Low, NA" | mcaide_apoe2 == "Mid, NA" | 
        mcaide_apoe2 == "High, NA" ~ NA_character_,
      mcaide_apoe2 == "Low, e2+" ~ "Low, e2+",
      mcaide_apoe2 == "Low, e3/e3" ~ "Low, e3/e3",
      mcaide_apoe2 == "Low, e4+" ~ "Low, e4+",
      mcaide_apoe2 == "Mid, e2+" ~ "Mid, e2+",
      mcaide_apoe2 == "Mid, e3/e3" ~ "Mid, e3/e3",
      mcaide_apoe2 == "Mid, e4+" ~ "Mid, e4+",
      mcaide_apoe2 == "High, e2+" ~ "High, e2+",
      mcaide_apoe2 == "High, e3/e3" ~ "High, e3/e3",
      mcaide_apoe2 == "High, e4+" ~ "High, e4+"
    ),
    mcaide_apoe2 = fct_relevel(mcaideXapoe2, "Mid, e3/e3", "Low, e2+", "Low, e3/e3",
                              "Low, e4+",
                              "Mid, e2+", "Mid, e4+", "High, e2+", "High, e3/e3",
                              "High, e4+")
  ) %>%
  mutate_at(
    vars(caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp, 
         caide_chol2, caide_chol3,     
      mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp, 
         mcaide_chol2, mcaide_chol3), as_factor
  ) %>% 
  # check for duplicates
  group_by(rid) %>% 
  arrange(rid, examdate) %>%
  slice_min(examdate) %>%
  ungroup() %>%
  dplyr::select(-c('examdate'))

janitor::get_dupes(caide_calc) # no dupes

caide_calc %>% count(is.na(apoe)) # NAs 195

write_csv(caide_calc, 
          "~/Desktop/ADNI/adni_caide_calc.csv")

str(caide_calc) # 2369 rows
# rid, dx_bl, cdrsb, i_age, i_education, i_bmi, i_systbp,
# rct20_num, rct20_no_imp, hypchol3, hypchol1, hypchol2

# caide_chol3, caide_chol1, caide_chol2, mcaide_chol3, mcaide_chol2, 

# apoe, gender, race 
# caide_age, caide_educ, caide_sex, caide_obesity, caide_sbp
# mcaide_age, mcaide_educ, mcaide_sex, mcaide_obesity, mcaide_sbp

# caide3, z_caide3, caide_apoe3, mcaide3, z_mcaide3, mcaide_apoe3
# caide1, z_caide1, caide_apoe1, mcaide1, z_mcaide1, mcaide_apoe1 
# caide2, z_caide2, caide_apoe2, mcaide2, z_mcaide2, mcaide_apoe2

# caide_missing3, caide_cat3, caideXapoe3, caide_missing1, caide_cat1, 
# caideXapoe1, mcaide_missing3, mcaide_cat3, mcaideXapoe3, mcaide_missing1, 
# mcaide_cat1, mcaideXapoe1, caide_missing2, caide_cat2, caideXapoe2, 
# mcaide_missing2, mcaide_cat2, mcaideXapoe2


##################################################################
# Make tables for caide and mcaide (without chol)
# Tables are stratified by RACE
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







