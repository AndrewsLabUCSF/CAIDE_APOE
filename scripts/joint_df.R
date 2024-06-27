library(tidyverse)
library(janitor)

## Import datasets 
nacc.raw <- read_csv('data/nacc.csv')
adni.raw <- read_tsv('data/adni.csv')


## ADNI 
#adni.raw %>% count(adrd)

adni <- adni.raw %>%
  select(ptid, origprot, age, ptgender, race, apoe, apoe_geno, pteducat, bmi, vsbpsys, total_c, htn, hld,
         i_bmi, i_vsbpsys, i_pteducat, i_total_c, i_htn, i_hld,
         i_caide_chol, starts_with('caide'), i_mcaide_chol, starts_with('mcaide'), 
         cdrsb, dx, naccudsd, naccetpr) %>%
  mutate(
    cohort = "ADNI", 
  ) %>%
  rename(
    hypchol = hld, 
    hyperten = htn, 
    gender = ptgender, 
    educ = pteducat, 
    bpsys = vsbpsys, 
    i_educ = i_pteducat, 
    i_bpsys = i_vsbpsys, 
    i_hyperten = i_htn, 
    i_hypchol = i_hld, 
  )
  
## NACC 
nacc <- nacc.raw %>% 
  select(NACCID, NACCAGE, SEX, race, apoe, apoe_geno, EDUC,  NACCBMI, BPSYS, hypchol, hyperten, 
         i_NACCBMI, i_BPSYS, i_EDUC, i_hypchol, i_hyperten, 
         starts_with('caide'), starts_with('mcaide'), 
         CDRSUM, dx, NACCUDSD, NACCETPR
         ) %>%
  select(-mcaide_cat, -mcaide_apoe, -caide_cat, -caide_apoe) %>%
  mutate(cohort = "NACC", i_caide_chol = caide_chol, i_mcaide_chol = caide_chol) %>%
  rename(
    ptid = NACCID, 
    age = NACCAGE, 
    gender = SEX, 
    educ = EDUC, 
    bmi = NACCBMI, 
    bpsys = BPSYS, 
    i_bmi = i_NACCBMI, 
    i_bpsys = i_BPSYS, 
    i_educ = i_EDUC, 
    cdrsb = CDRSUM, 
    naccudsd = NACCUDSD, 
    naccetpr = NACCETPR
    ) %>%
  mutate(
    gender = ifelse(gender == 1, "Male", "Female")
  )  

setdiff(names(adni), names(nacc))
get_na_by_cols(adni)

## Joint dataset
joint <- bind_rows(
  adni, nacc
) %>%
  mutate(
    race = ifelse(race == "Non-Hispanic White", "NHW", race)
  )  %>%
  mutate(
    ## CAIDE
    z_caide = scale(caide)[,1],
    z_caide_nosex = scale(caide_nosex)[,1],
    caide_cat = case_when(
      between(caide, 0, (mean(caide, na.rm = T) - sd(caide, na.rm = T))) ~ 'low',
      between(caide, (mean(caide, na.rm = T) - sd(caide, na.rm = T)), (mean(caide, na.rm = T) + sd(caide, na.rm = T))) ~ 'mid', 
      between(caide, (mean(caide, na.rm = T) + sd(caide, na.rm = T)), 14) ~ 'high',
      TRUE ~ NA_character_
    ), 
    caide_cat = fct_relevel(caide_cat, 'low', "mid",  'high'),
    caide_apoe = glue("{caide_cat}_{apoe}"), 
    caide_apoe = fct_relevel(caide_apoe, "mid_e3/e3", "low_e2+", "low_e3/e3", "low_e4+", 
                             "mid_e2+", "mid_e4+", "high_e2+", "high_e3/e3", "high_e4+")
  ) %>%
  mutate(
    ## mCAIDE
    z_mcaide = scale(mcaide)[,1],
    z_mcaide_nosex = scale(mcaide_nosex)[,1],
    mcaide_cat = case_when(
      between(mcaide, 0, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T))) ~ 'low',
      between(mcaide, (mean(mcaide, na.rm = T) - sd(mcaide, na.rm = T)), (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T))) ~ 'mid', 
      between(mcaide, (mean(mcaide, na.rm = T) + sd(mcaide, na.rm = T)), 14) ~ 'high',
      TRUE ~ NA_character_
    ), 
    mcaide_cat = fct_relevel(mcaide_cat, 'low', "mid",  'high'),
    mcaide_apoe = glue("{mcaide_cat}_{apoe}"), 
    mcaide_apoe = fct_relevel(mcaide_apoe, "mid_e3/e3", "low_e2+", "low_e3/e3", "low_e4+", 
                              "mid_e2+", "mid_e4+", "high_e2+", "high_e3/e3", "high_e4+")
  ) 

## Export 
write_tsv(joint, "data/joint_df.tsv.gz")

## Table
tab_rf <- joint %>% 
  dplyr::select(race, cohort, gender, age, i_educ, i_bmi, i_hyperten, i_hypchol, 
                caide, caide_cat, mcaide, mcaide_cat, apoe, apoe_geno, naccudsd, naccetpr) %>%
  mutate(
    naccudsd = fct_recode(as.factor(naccudsd), 'CU' = '1', "MCI" = "3", "ADRD" = "4"),
    naccetpr = fct_recode(as.factor(naccetpr), 
                          'CU' = '88', "AD" = "1", "LBD" = "2", "PSP" = '4',  
                          'CBD' = "5", "FTLD" = '6', "FTLD" = '7', 'VCID' = '8'),
    naccetpr = fct_relevel(naccetpr),
    # SEX = as.factor(SEX),
    # SEX = fct_recode(SEX,  'Male' = '1', "Female" = "2")
  ) %>%
  tbl_summary(by = race,
              statistic = list(all_continuous() ~ "{mean} ({sd})"), 
              label = list(gender ~ "Gender",
                           age ~ "Age",
                           i_educ ~ "Education",
                           race ~ "Race",
                           naccudsd ~ "Diagnosis",
                           naccetpr ~ "ADRD",
                           i_bmi ~ "BMI",
                           i_hyperten ~ "Hypertension",
                           i_hypchol ~ "Hypercholesterolemia",
                           caide ~ "CAIDE",
                           mcaide ~ "mCAIDE",
                           apoe ~ "APOE", 
                           cohort ~ "Cohort", 
                           caide_cat ~ "CAIDE",
                           mcaide_cat ~ "mCAIDE"
                           
              )
              ) %>%
  add_n() 

tab_rf










































