library(tidyverse)
library(janitor)

## Import datasets 
nacc.raw <- read_csv('data/nacc.csv')
adni.raw <- read_tsv('data/adni.csv')
adprs.raw <- read_tsv('~/Downloads/results_newfam/out/score/out_pgs.txt.gz', guess_max = 60000) %>%
  janitor::clean_names()
pcs <- read_tsv('~/Downloads/results_newfam/out/score/out_popsimilarity.txt.gz', guess_max = 60000)  %>%
  janitor::clean_names()
test <- read_table("~/Downloads/adgc_META_pst_eff_a1_b0.5_phiauto_chrAll_noAPOE_scores.txt", guess_max = 60000)  %>%
  janitor::clean_names()

## Genetics
prs <- left_join(adprs.raw, pcs, by = c('sampleset', 'iid')) %>%
  filter(., str_detect(iid, "NACC")) %>%
  mutate(
    iid = str_replace(iid, "^\\d+_", ""),
    iid = str_replace(iid, "^[A-Za-z0-9]+_", "")
  ) %>%
  select(-sampleset, -pgs) 

filter(test, str_detect(iid, "NACC"))
filter(prs, str_detect(iid, "NACC"))

## ADNI 
adni <- adni.raw %>%
  select(ptid, origprot, age, ptgender, race, apoe, apoe_geno, pteducat, bmi, 
         vsbpsys, total_c, htn, hld,
         i_pteducat, i_bmi, i_vsbpsys, i_total_c, i_htn, i_hld, i_caide_chol,
         starts_with('caide'), i_mcaide_chol, starts_with('mcaide'), 
         starts_with("m2caide"),
         cdrsb, dx, naccudsd, naccetpr,
         #add in unimputed, dichotimize the adni cholesterol variable
         caide_chol, mcaide_chol
         ) %>%
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
  select(NACCID, NACCAGE, SEX, race, apoe, apoe_geno,   
         EDUC, NACCBMI, BPSYS, hypchol, hyperten, 
         i_EDUC, i_NACCBMI, i_BPSYS,  i_hypchol, i_hyperten, 
         starts_with('caide'), starts_with('mcaide'), 
         starts_with('m2caide'),
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

## Joint dataset
joint <- bind_rows(
  adni, nacc
) %>%
  left_join(prs, by = c("ptid" = "iid")) %>%
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
                             "mid_e2+", "mid_e4+", "high_e2+", "high_e3/e3", "high_e4+"),
    
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
  ) %>%
  group_by(most_similar_pop) %>%
  mutate(
    ## PRS 
    z_prs = scale(z_norm2)[,1], 
    prs_cat = case_when(
      between(z_prs, -Inf, (mean(z_prs, na.rm = T) - sd(z_prs, na.rm = T))) ~ 'low',
      between(z_prs, (mean(z_prs, na.rm = T) - sd(z_prs, na.rm = T)), (mean(z_prs, na.rm = T) + sd(z_prs, na.rm = T))) ~ 'mid', 
      between(z_prs, (mean(z_prs, na.rm = T) + sd(z_prs, na.rm = T)), Inf) ~ 'high',
      TRUE ~ NA_character_
    ), 
    prs_cat = fct_relevel(prs_cat, 'low', "mid",  'high'),
    mcaide_prs = glue("{mcaide_cat}_{prs_cat}"), 
    mcaide_prs = ifelse(str_detect(mcaide_prs, "NA"), NA_character_, mcaide_prs),
    mcaide_prs = fct_relevel(mcaide_prs, "mid_mid", "low_low", "low_mid", "low_high", 
                              "mid_low", "mid_high", "high_low", "high_mid", "high_high")
    
  ) %>%
  ungroup() %>% 
  mutate(
    #m2CAIDE
    #create zscore for caide
    z_m2caide = scale(m2caide)[,1],
    
    m2caide_cat = case_when(
      between(m2caide, 0, (mean(m2caide, na.rm = T) - sd(m2caide, na.rm = T))) ~ 'low',
      between(m2caide, (mean(m2caide, na.rm = T) - sd(m2caide, na.rm = T)), (mean(m2caide, na.rm = T) + sd(m2caide, na.rm = T))) ~ 'mid', 
      between(m2caide, (mean(m2caide, na.rm = T) + sd(m2caide, na.rm = T)), 14) ~ 'high',
      TRUE ~ NA_character_
    ),
    m2caide_cat = fct_relevel(m2caide_cat, 'low', "mid",  'high'),
    m2caide_apoe = glue("{m2caide_cat}_{apoe}"), 
    m2caide_apoe = fct_relevel(m2caide_apoe, "mid_e3/e3", "low_e2+", "low_e3/e3", "low_e4+", 
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

raw_prs <- bind_rows(
  adni, nacc
) %>%
  left_join(filter(test, str_detect(iid, "NACC")), by = c("ptid" = "iid")) %>%
  select(ptid, origprot, age, race, scoresum)

summary(raw_prs)

norm_prs <- bind_rows(
  adni, nacc
) %>%
  left_join(prs, by = c("ptid" = "iid")) %>%
  select(ptid, origprot, age, race, z_norm2, most_similar_pop)

summary(norm_prs)

left_join(raw_prs, norm_prs) %>%
  filter(!is.na(scoresum), is.na(z_norm2))

filter(prs, str_detect(iid, "NACC002431"))
filter(adprs.raw, str_detect(iid, "NACC")) %>%
  separate(iid, c('v1', 'v2', 'v3')) %>%
  mutate(first_char = str_sub(v2, 1, 4)) %>%
  count(first_char)
  
filter(adprs.raw, str_detect(iid, "NACC")) %>%
  mutate(
    iid = str_replace(iid, "^\\d+_", ""),
    iid = str_replace(iid, "^[A-Za-z0-9]+_", "")
  ) %>%
  mutate(first_char = str_sub(iid, 1, 2)) %>%
  filter(first_char == "AD") %>%
  count(first_char) %>% print(n = Inf)

prs <- left_join(adprs.raw, pcs, by = c('sampleset', 'iid')) %>%
  mutate(
    iid = str_replace(iid, "0_", "")
  ) %>%
  select(-sampleset, -pgs) 

  
filter(nacc, str_detect(ptid, "09AD14012_NACC929764"))
filter(adprs.raw, str_detect(iid, "08AD11427"))


nacc %>%
  mutate(first_char = str_sub(ptid, 1, 4)) 
  































