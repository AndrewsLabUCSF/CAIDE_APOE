library(tidyverse)
library(janitor)
# library(gtsummary)
library(missForest)

`%nin%` = negate(`%in%`)

setwd("~/gitcode/CAIDE_APOE")

# NACC
## Import datasets
adrc.raw <- read_csv("resources/NACC/ADRC_all_9_16_2022.csv")
adgc.raw <- read_tsv('resources/NACC/ADGC_NACCIDs.txt', col_names = F) %>%
  rename(NACCID = X1) %>%
  mutate(adgc = T)
nacc.raw <- read_csv('resources/NACC/investigator_nacc58.csv')

## Extract APOE genotypes
adrc <- adrc.raw %>%
  unite(apoe_geno, APOE1, APOE2, sep = '/') %>%
  mutate(apoe_geno = ifelse(apoe_geno == "NA/NA", NA, apoe_geno), 
         gwas = ifelse(is.na(GWAS), FALSE, TRUE), 
         wgs = ifelse(is.na(WGS), FALSE, TRUE)) %>%
  dplyr::select(NACCID, apoe_geno, gwas, wgs) %>%
  distinct() %>%
  mutate(
    apoe_geno = fct_recode(apoe_geno,
                            '22' = "E2/E2", '23' = "E2/E3", '23' = "E3/E2", 
                            '33' = "E3/E3", 
                            '24' = "E2/E4", '24' = "E4/E2", 
                            '34' = "E4/E3", '34' = "E3/E4", '44' = "E4/E4"
                            ),
    apoe = case_when(
      apoe_geno %in% c("22", "23") ~ "e2+",
      apoe_geno %in% c("33") ~ "e3/e3",
      apoe_geno %in% c("34", "24", "44") ~ "e4+",
      TRUE ~ NA_character_
    ), 
    apoe = fct_relevel(apoe, "e3/e3"),
    apoe4 = case_when(
      apoe_geno %in% c("22", "23", "33") ~ "e4-",
      apoe_geno %in% c("24", "34", "44") ~ "e4+",
      TRUE ~ NA_character_
    ), 
    apoe4 = fct_relevel(apoe4, "e4-"),
  )

## NACC variables of interest
nacc_df <- select(nacc.raw, NACCID, NACCADC, VISITYR, PACKET, NACCVNUM, NACCAVST,
       BIRTHYR, SEX, EDUC, NACCAGE, NACCAGEB, NACCDIED, NACCYOD, 
       RACE, HISPANIC, NACCNIHR, 
       HEIGHT, WEIGHT, NACCBMI, 
       DIABETES, DIABTYPE, DIABET, NACCDBMD,
       BPSYS, BPDIAS, HYPERTEN, HYPERT, NACCAHTN,
       HYPERCHO, HYPCHOL, NACCLIPL, 
       CBSTROKE, STROKMUL, HXSTROKE, PREVSTK, STKIMAG, 
       CVAFIB, AFIBRILL, 
       TOBAC30, TOBAC100, SMOKYRS, PACKSPER, QUITSMOK,
       NACCTBI, 
       INSOMN, HYPOSOM, 
       NACCGDS, 
       NPSYLAN, 
       NACCMMSE, CDRSUM, LOGIMEM,
       NORMCOG, NACCTMCI, DEMENTED, NACCUDSD, # Prevelant
       NACCNORM, NACCMCII, NACCIDEM, # Incident
       DECAGE, NACCETPR, 
       NACCADMU, NACCFTDM
       ) 

## clean data 
### append genotypes, specificing missing values
nacc_clean <- nacc_df %>%
  left_join(adrc, by = 'NACCID') %>%
  left_join(adgc.raw, by = 'NACCID') %>%
  mutate(
    ## Missing data 
    NACCBMI = na_if(NACCBMI, 888.8),
    NACCBMI = na_if(NACCBMI, -4), 
    NACCTBI = na_if(NACCTBI, 9), 
    NACCTBI = na_if(NACCTBI, -4), 
    BPSYS = na_if(BPSYS, -4), 
    BPSYS = na_if(BPSYS, 777), 
    BPSYS = na_if(BPSYS, 888),
    BPDIAS = na_if(BPDIAS, -4), 
    BPDIAS = na_if(BPDIAS, 777), 
    BPDIAS = na_if(BPDIAS, 888), 
    NACCMMSE = ifelse(between(NACCMMSE, 0, 30), NACCMMSE, NA), 
    EDUC = na_if(EDUC, 99), 
    HEIGHT = na_if(HEIGHT, -4), 
    HEIGHT = na_if(HEIGHT, 88.8), 
    WEIGHT = na_if(WEIGHT, -4), 
    WEIGHT = na_if(WEIGHT, 888.0), 
    DIABETES = na_if(DIABETES, -4), 
    DIABETES = na_if(DIABETES, 9), 
    DIABTYPE = na_if(DIABTYPE, -4), 
    DIABTYPE = na_if(DIABTYPE, 9), 
    DIABTYPE = na_if(DIABTYPE, 8), 
    DIABET = na_if(DIABET, 9), 
    DIABET = na_if(DIABET, -4), 
    NACCDBMD = na_if(NACCDBMD, -4), 
    NACCAHTN = na_if(NACCAHTN, -4), 
    HYPERTEN = na_if(HYPERTEN, -4),
    HYPERTEN = na_if(HYPERTEN, 9),
    HYPERT = na_if(HYPERT, -4),
    HYPERT = na_if(HYPERT, 8),
   HYPERCHO = na_if(HYPERCHO, -4),
   HYPERCHO = na_if(HYPERCHO, 9),
   HYPCHOL = na_if(HYPCHOL, -4),
   HYPCHOL = na_if(HYPCHOL, 8),
   NACCLIPL = na_if(NACCLIPL, -4),
   CBSTROKE = na_if(CBSTROKE, -4),
   CBSTROKE = na_if(CBSTROKE, 9),
   STROKMUL = na_if(STROKMUL, 9),
   STROKMUL = na_if(STROKMUL, -4),
   STROKMUL = na_if(STROKMUL, 8),
   HXSTROKE = na_if(HXSTROKE, -4),
   PREVSTK = na_if(PREVSTK, -4),
   STKIMAG = na_if(STKIMAG, 8),
   STKIMAG = na_if(STKIMAG, 9),
   STKIMAG = na_if(STKIMAG, -4),
   CVAFIB = na_if(CVAFIB, -4),
   CVAFIB = na_if(CVAFIB, 9),
   AFIBRILL = na_if(AFIBRILL, -4),
   AFIBRILL = na_if(AFIBRILL, 8),
   TOBAC30 = na_if(TOBAC30, -4),
   TOBAC30 = na_if(TOBAC30, 9),
   TOBAC100 = na_if(TOBAC100, -4),
   TOBAC100 = na_if(TOBAC100, 9),
   SMOKYRS = ifelse(between(SMOKYRS, 0, 87), SMOKYRS, NA),
   PACKSPER = na_if(PACKSPER, -4),
   PACKSPER = na_if(PACKSPER, 8),
   PACKSPER = na_if(PACKSPER, 9),
   QUITSMOK = ifelse(between(QUITSMOK, 0, 110), QUITSMOK, NA),
   INSOMN = na_if(INSOMN, -4),
   INSOMN = na_if(INSOMN, 9),
   HYPOSOM = na_if(HYPOSOM, -4),
   HYPOSOM = na_if(HYPOSOM, 8),
   NACCGDS = ifelse(between(NACCGDS, 0, 15), NACCGDS, NA),
   NPSYLAN = na_if(NPSYLAN, -4),
   DECAGE = ifelse(between(DECAGE, 0, 110), DECAGE, NA),
   NACCYOD = na_if(NACCYOD, 8888),
  ) %>%
  mutate(
    race = case_when(
      NACCNIHR == 1 & HISPANIC == 0 ~ "NHW",
      NACCNIHR == 1 & HISPANIC == 1 ~ "Hispanic",
      NACCNIHR == 2 ~ "Black",
      NACCNIHR == 5 ~ "Asian",
      TRUE ~ NA_character_
    ),
    race = fct_relevel(race, 'NHW')
  ) 

nacc_wrangle <- nacc_clean %>%
  filter(NACCVNUM == 1) %>%
  filter(NACCAGE >= 55) %>%
  filter(is.na(DECAGE) | DECAGE >= 55)  %>%
  filter(NACCNIHR %in% c(1,2,5)) %>%
  filter(NACCADMU == 0) %>%
  filter(NACCFTDM == 0) %>%
  filter(NACCETPR %in% c(1, 2, 4, 5, 6, 7, 8, 88)) %>%
  filter(NACCUDSD %nin% c(2)) %>%
  filter(!is.na(race)) %>%
  filter(!is.na(apoe)) %>%
  mutate(
    hypchol = case_when(
      HYPCHOL == 0 ~ 0, 
      HYPCHOL == 1 ~ 1, 
      is.na(HYPCHOL) & HYPERCHO == 0 ~ 0, 
      is.na(HYPCHOL) & HYPERCHO == 1 ~ 1, 
      TRUE ~ NA_real_
    ), 
    diabetes = case_when(
      DIABET == 0 ~ 0, 
      DIABET %in% c(1, 2, 3) ~ 1,
      is.na(DIABET) & DIABETES == 0 ~ 0, 
      is.na(DIABET) & DIABETES %in% c(1,2) ~ 1, 
      TRUE ~ NA_real_
    ), 
    stroke = case_when(
      PREVSTK == 0 ~ 0, 
      PREVSTK == 1 ~ 1,
      is.na(PREVSTK) & CBSTROKE == 0 ~ 0, 
      is.na(PREVSTK) & CBSTROKE %in% c(1,2) ~ 1,
      TRUE ~ NA_real_
    ), 
    hyperten = case_when(
      HYPERT == 0 ~ 0, 
      HYPERT == 1 ~ 1, 
      is.na(HYPERT) & HYPERTEN == 0 ~ 0, 
      is.na(HYPERT) & HYPERTEN %in% c(1,2) ~ 1,
      TRUE ~ NA_real_
    ), 
    insomnia = case_when(
      HYPOSOM == 0 ~ 0, 
      HYPOSOM == 1 ~ 1,
      is.na(HYPOSOM) & INSOMN == 0 ~ 0, 
      is.na(HYPOSOM) & INSOMN %in% c(1,2) ~ 1,
      TRUE ~ NA_real_
    ),
    afib = case_when(
      AFIBRILL == 0 ~ 0, 
      AFIBRILL == 1 ~ 1,
      is.na(AFIBRILL) & CVAFIB == 0 ~ 0, 
      is.na(AFIBRILL) & CVAFIB %in% c(1,2) ~ 1,
      TRUE ~ NA_real_
    ),
    smk = case_when(
      TOBAC100 == 0 ~ 0, 
      TOBAC100 == 0 & TOBAC30 == 1 ~ 0, 
      TOBAC100 == 1 & TOBAC30 == 0 ~ 0, 
      TOBAC100 == 1 & TOBAC30 == 1 ~ 1, 
      TRUE ~ NA_real_
    ), 
    dx = ifelse(NACCUDSD %in% c(3,4), 1, 0),
    incident_ci = case_when(
      NACCIDEM == 1 ~ 'AD',
      NACCNORM == 1 ~ 'CN', 
      NACCMCII == 1 ~ 'MCI', 
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate_at(
    vars(SEX, DIABETES, DIABTYPE, DIABET, NACCDBMD, NACCAHTN, HYPERTEN, HYPERT, 
         HYPERCHO, HYPCHOL, NACCLIPL, CBSTROKE, STROKMUL, HXSTROKE, PREVSTK, STKIMAG, 
         CVAFIB, AFIBRILL, TOBAC30, TOBAC100, PACKSPER, NACCTBI, INSOMN, HYPOSOM, NPSYLAN, 
         apoe_geno, hypchol, diabetes, stroke, hyperten, insomnia, afib, smk), as_factor
  )
  

### impute missing data 
library(doParallel)
registerDoParallel(cores=3) 
getDoParWorkers()
require(doRNG)
registerDoRNG(seed = 333)

nacc_ximp <- nacc_wrangle %>% 
  select(NACCAGE, race, CDRSUM, apoe_geno, SEX, HEIGHT, WEIGHT, NACCBMI, BPSYS, NACCTBI, NACCGDS, 
         BPDIAS, EDUC, hypchol, diabetes, stroke, hyperten, insomnia, afib, smk) %>%
  as.data.frame() %>%
  missForest(xmis = ., variablewise = TRUE, verbose = TRUE, 
             xtrue = filter(., complete.cases(.)), parallelize = 'forests')

nacc_imp <- nacc_ximp %>%
  magrittr::extract2(1) %>%
  select(-NACCAGE, -race, -CDRSUM, -apoe_geno, -SEX) %>%
  magrittr::set_colnames(paste0("i_", colnames(.))) %>%
  as_tibble() %>%
  mutate(
    i_EDUC = round(i_EDUC), 
    i_NACCBMI = round(i_NACCBMI, 1)
  )

nacc <- nacc_wrangle %>%
  bind_cols(nacc_imp)

## NACC Last visit 
## NACC variables of interest
nacc_lv2 <- nacc_clean %>%
  filter(NACCID %in% nacc$NACCID) %>%
  select(., NACCID, NACCADC, VISITYR, PACKET, NACCVNUM,
                      NACCAGE,
                      NACCMMSE, CDRSUM,
                      NACCUDSD, # Prevelant
                      DECAGE, NACCETPR
         ) %>%
  ## Retain only information from last visit
  filter(NACCETPR %nin% c(99)) %>%
  group_by(NACCID) %>%
  filter(NACCVNUM == which.max(NACCVNUM)) %>%
  ungroup() %>%
  ## Age of Onset, for diagnosis, age at which symptoms began, for controls, age at last examination.
  filter(NACCUDSD != 2) %>%
  mutate(
    aoo = case_when(
      NACCETPR %in% c(1:30) ~ DECAGE,
      NACCETPR == 88 ~ NACCAGE,
      TRUE ~ NA_real_
    ),
    dx = ifelse(NACCUDSD %in% c(3,4), 1, 0),
    ad_age = case_when(
      dx == 0 ~ log(1-((NACCAGE-54.5)/(110.5-54.5))) - 0.5,
      dx == 1 ~ -log((NACCAGE-54.5)/(110.5-54.5)) + 0.5
    ), 
    sdx = case_when(
      dx == 0 ~ log(1-((aoo-54.5)/(110.5-54.5))) - 0.5,
      dx == 1 ~ -log((aoo-54.5)/(110.5-54.5)) + 0.5
    )
  ) %>%
  magrittr::set_colnames(paste0(colnames(.), "_lv")) %>%
  rename(NACCID = NACCID_lv)

nacc_lv_wrangle <-  nacc_clean %>%
  filter(NACCID %in% nacc$NACCID) %>%
  select(., NACCID, race, NACCAVST, NACCADC, VISITYR, PACKET, NACCVNUM, NACCAVST, NACCDIED,
         NACCAGE, 
         NACCMMSE, CDRSUM,
         NACCUDSD, # Prevelant,
         NACCIDEM, NACCMCII,
         DECAGE, NACCETPR
  ) %>%
  ## Retain only information from last visit
  filter(NACCETPR %nin% c(99)) %>%
  filter(NACCUDSD != 2) %>%
  mutate(
    ADRD = case_when(
      NACCETPR == 88 ~ 0,
      NACCETPR %in% c(1, 2, 6, 7, 8) ~ 1, 
      NACCETPR %in% c(3, 4, 5, 9:30) ~ 88,
      TRUE ~ NA_real_
    )
  )

## CI at baseline or CI not due ADRD 
exclude_ID <- filter(nacc_lv_wrangle, (NACCVNUM == 1 & NACCUDSD %in% c(3, 4)) | ADRD == 88) %>% distinct(NACCID) %>% pull(NACCID)
exclude_ADRD_ID <- filter(nacc_lv_wrangle, ADRD %in% c(1, 88)) %>% distinct(NACCID) %>% pull(NACCID)

adrd_aoo <- nacc_lv_wrangle %>% 
  filter(NACCID %nin% exclude_ID) %>%
  group_by(NACCID, NACCUDSD) %>%
  summarise(
    NACCVNUM = min(NACCVNUM),
    NACCAGE = min(NACCAGE), 
    # race = first(race),
  ) %>%
  filter(NACCUDSD != 1) %>% 
  ungroup() %>%
  mutate(
    NACCUDSD = as.character(NACCUDSD),
    NACCUDSD = fct_recode(NACCUDSD, "MCI" = "3", "Dementia" = "4")
  ) %>%
  select(-NACCVNUM) %>%
  pivot_wider(names_from = NACCUDSD, values_from = NACCAGE) %>%
  mutate(
    dx_aoo = 1, 
  ) %>%
  rename(mci_aoo = MCI, dementia_aoo = Dementia)

ctrl_aoo <- nacc_lv_wrangle %>% 
  filter(NACCID %nin% exclude_ADRD_ID) %>% 
  group_by(NACCID) %>%
  summarise(
    NACCVNUM = max(NACCVNUM),
    NACCAGE = max(NACCAGE), 
    # race = first(race),
  ) %>%
  ungroup() %>%
  # filter(NACCVNUM != 1) %>%
  select(-NACCVNUM) %>%
  mutate(
    dx_aoo = 0
  ) %>%
  rename(ctrl_aoo = NACCAGE)

nacc_lv <- bind_rows(
  adrd_aoo, ctrl_aoo
) %>%
  mutate(
    aoo = case_when(
      dx_aoo == 1 & is.na(dementia_aoo) ~ mci_aoo,
      dx_aoo == 1 & mci_aoo <= dementia_aoo ~ mci_aoo, 
      dx_aoo == 1 & is.na(mci_aoo) ~ dementia_aoo,
      dx_aoo == 1 & dementia_aoo < mci_aoo ~ dementia_aoo, 
      dx_aoo == 0 ~ ctrl_aoo
    )
  )


## Dementia risk scores
### CAIDE
caide <- nacc %>%
  mutate(
    caide_age = case_when(
      NACCAGE <= 46 ~ 0, 
      NACCAGE >= 47 & NACCAGE < 53 ~ 3, 
      NACCAGE >= 53 ~ 4, 
      TRUE ~ NA_real_
    ), 
    caide_educ = case_when(
      i_EDUC <= 6 ~ 3, 
      i_EDUC >= 7 & i_EDUC <= 9 ~ 2, 
      i_EDUC > 9 ~ 0, 
      TRUE ~ NA_real_
    ), 
    caide_sex = case_when(
      SEX == 1 ~ 1, 
      SEX == 2 ~ 0, 
      TRUE ~ NA_real_
    ), 
    caide_bmi = case_when(
      i_NACCBMI < 30 ~ 0, 
      i_NACCBMI >= 30 ~ 2, 
      TRUE ~ NA_real_
    ),
    caide_sbp = case_when(
      i_BPSYS < 140 ~ 0, 
      i_BPSYS >= 140 ~ 2, 
      TRUE ~ NA_real_
    ),
    caide_chol = case_when(
      i_hypchol == 0 ~ 0, 
      i_hypchol == 1 ~ 2, 
      TRUE ~ NA_real_
    ), 
  ) %>%
  rowwise() %>%
  mutate(
    caide = sum(caide_age, caide_educ, caide_sex, caide_bmi, caide_sbp, caide_chol, na.rm = F), 
    caide_nosex = sum(caide_age, caide_educ, caide_bmi, caide_sbp, caide_chol, na.rm = F), 
    caide_missing = sum(is.na(NACCAGE), is.na(EDUC), is.na(SEX), 
                        is.na(NACCBMI), is.na(hypchol), is.na(BPSYS)),
  ) %>%
  ungroup() %>%
  mutate(
    z_caide = scale(caide)[,1],
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
  mutate_at(
    vars(caide_age, caide_educ, caide_sex, caide_bmi, caide_sbp, caide_chol), as_factor
  ) %>%
  select(NACCID, starts_with('caide'), z_caide, caide_cat)
  
### mCAIDE
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8483620/
mcaide <- nacc %>%
  mutate(
    mcaide_age = case_when(
      NACCAGE <= 64 ~ 0, 
      NACCAGE >= 65 & NACCAGE < 73 ~ 1, 
      NACCAGE >= 73 ~ 2, 
      TRUE ~ NA_real_
    ), 
    mcaide_educ = case_when(
      i_EDUC <= 11 ~ 2, 
      i_EDUC >= 12 & i_EDUC <= 16 ~ 1, 
      i_EDUC > 16 ~ 0, 
      TRUE ~ NA_real_
    ), 
    mcaide_sex = case_when(
      SEX == 1 ~ 1, 
      SEX == 2 ~ 0, 
      TRUE ~ NA_real_
    ), 
    mcaide_bmi = case_when(
      i_NACCBMI <= 30 ~ 0, 
      i_NACCBMI > 30 ~ 2, 
      TRUE ~ NA_real_
    ),
    mcaide_sbp = case_when(
      i_BPSYS < 140 ~ 0, 
      i_BPSYS >= 140 ~ 2, 
      TRUE ~ NA_real_
    ),
    mcaide_chol = case_when(
      i_hypchol == 0 ~ 0, 
      i_hypchol == 1 ~ 2, 
      TRUE ~ NA_real_
    ), 
  ) %>%
  rowwise() %>%
  mutate(#insert new mcaide without age nor sex here. 
    mcaide = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_bmi, mcaide_sbp, mcaide_chol, na.rm = F), 
    mcaide_nosex = sum(mcaide_age, mcaide_educ, mcaide_bmi, mcaide_sbp, mcaide_chol, na.rm = F), 
    mcaide_missing = sum(is.na(NACCAGE), is.na(EDUC), is.na(SEX), 
                        is.na(NACCBMI), is.na(hypchol), is.na(BPSYS)),
  ) %>%
  ungroup() %>%
  mutate(
    z_mcaide = scale(mcaide)[,1],
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
  mutate_at(
    vars(mcaide_age, mcaide_educ, mcaide_sex, mcaide_bmi, mcaide_sbp, mcaide_chol), as_factor
  ) %>%
  select(NACCID, starts_with('mcaide'), z_mcaide, mcaide_cat)

## CogDrisk Score 
cogdrisk <- nacc %>%
  mutate(
    cogdrisk_age = case_when(
      SEX == 1 & between(NACCAGE, 0, 64) ~ 0, 
      SEX == 1 & between(NACCAGE, 65, 69) ~ 6, 
      SEX == 1 & between(NACCAGE, 70, 74) ~ 8, 
      SEX == 1 & between(NACCAGE, 75, 79) ~ 13, 
      SEX == 1 & between(NACCAGE, 80, 84) ~ 17, 
      SEX == 1 & between(NACCAGE, 85, 89) ~ 20, 
      SEX == 1 & between(NACCAGE, 90, Inf) ~ 22, 
      SEX == 2 & between(NACCAGE, 0, 64) ~ 0, 
      SEX == 2 & between(NACCAGE, 65, 69) ~ 4, 
      SEX == 2 & between(NACCAGE, 70, 74) ~ 7, 
      SEX == 2 & between(NACCAGE, 75, 79) ~ 11, 
      SEX == 2 & between(NACCAGE, 80, 84) ~ 15, 
      SEX == 2 & between(NACCAGE, 85, 89) ~ 19, 
      SEX == 2 & between(NACCAGE, 90, Inf) ~ 23, 
      TRUE ~ NA_real_
    ), 
    cogdrisk_educ = case_when(
      i_EDUC < 8 ~ 4, 
      i_EDUC >= 8 & i_EDUC <= 11 ~ 2, 
      i_EDUC > 11 ~ 0, 
      TRUE ~ NA_real_
    ), 
    cogdrisk_bmi = case_when(
      between(i_NACCBMI, 10, 18.49) ~ 2,
      between(i_NACCBMI, 18.5, 24.99) ~ 0,
      between(i_NACCBMI, 25, 29.99) ~ 1,
      between(i_NACCBMI, 30, 100) ~ 3,
      TRUE ~ NA_real_
    ),
    cogdrisk_chol = case_when(
      i_hypchol == 0 ~ 0, 
      i_hypchol == 1 ~ 3, 
      TRUE ~ NA_real_
    ), 
    cogdrisk_diab = case_when(
      i_diabetes == 0 ~ 0, 
      SEX == 1 & i_diabetes == 1 ~ 2, 
      SEX == 2 & i_diabetes == 1 ~ 3, 
      TRUE ~ NA_real_
    ), 
    cogdrisk_stroke = case_when(
      i_stroke == 0 ~ 0, 
      i_stroke == 1 ~ 2,
      TRUE ~ NA_real_
    ), 
    cogdrisk_tbi = case_when(
      i_NACCTBI == 0 ~ 0, 
      i_NACCTBI == 1 ~ 2,
      TRUE ~ NA_real_
    ), 
    cogdrisk_bp = case_when(
      i_hyperten == 1 ~ 1, 
      i_BPSYS >= 130 ~ 1, 
      i_BPDIAS >= 80 ~ 1, 
      i_hyperten == 0 ~ 0, 
      i_BPSYS < 130 ~ 0, 
      i_BPDIAS < 80 ~ 0, 
      TRUE ~ NA_real_
    ), 
    # cogdrisk_insom = case_when(
    #   i_insomnia == 0 ~ 0, 
    #   i_insomnia == 1 ~ 2,
    #   TRUE ~ NA_real_
    # ),
    cogdrisk_dep = case_when(
      between(i_NACCGDS, 0, 9) ~ 0, 
      between(i_NACCGDS, 10, 15) ~ 3, 
      TRUE ~ NA_real_
    ), 
    cogdrisk_stroke = case_when(
      i_stroke == 0 ~ 0, 
      i_stroke == 1 ~ 2,
      TRUE ~ NA_real_
    ), 
    cogdrisk_afib = case_when(
      i_afib == 0 ~ 0, 
      i_afib == 1 ~ 2,
      TRUE ~ NA_real_
    ), 
    cogdrisk_smk = case_when(
      i_smk == 0 ~ 0, 
      i_smk == 1 ~ 1, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    cogdrisk = sum(cogdrisk_age, cogdrisk_educ, cogdrisk_bmi, cogdrisk_chol, 
                   cogdrisk_bmi, cogdrisk_chol, cogdrisk_diab, cogdrisk_stroke, 
                   cogdrisk_tbi, cogdrisk_bp, cogdrisk_dep, # cogdrisk_insom,
                   cogdrisk_smk), 
    cogdrisk_missing = sum(
      is.na(NACCAGE), is.na(EDUC), is.na(SEX), is.na(NACCBMI), is.na(hypchol), 
      is.na(diabetes), is.na(stroke), is.na(NACCTBI), is.na(NACCTBI), is.na(hyperten), 
      is.na(NACCGDS), is.na(stroke), is.na(afib), is.na(smk)
    )
  ) %>%
  ungroup() %>%
  mutate(
    z_cogdrisk = scale(cogdrisk)[,1],
    cogdrisk_cat = case_when(
      between(cogdrisk, 0, (mean(cogdrisk, na.rm = T) - sd(cogdrisk, na.rm = T))) ~ 'low',
      between(cogdrisk, (mean(cogdrisk, na.rm = T) - sd(cogdrisk, na.rm = T)), (mean(cogdrisk, na.rm = T) + sd(cogdrisk, na.rm = T))) ~ 'mid', 
      between(cogdrisk, (mean(cogdrisk, na.rm = T) + sd(cogdrisk, na.rm = T)), 50) ~ 'high',
      TRUE ~ NA_character_
    ), 
    cogdrisk_cat = fct_relevel(cogdrisk_cat, 'low', "mid",  'high'), 
    cogdrisk_apoe = glue("{cogdrisk_cat}_{apoe}"), 
    cogdrisk_apoe = fct_relevel(cogdrisk_apoe, "mid_e3/e3", "low_e2+", "low_e3/e3", "low_e4+", 
                              "mid_e2+", "mid_e4+", "high_e2+", "high_e3/e3", "high_e4+")
  ) %>%
  mutate_at(
    vars(cogdrisk_age, cogdrisk_educ, cogdrisk_bmi, cogdrisk_chol, cogdrisk_bp, 
         cogdrisk_diab, cogdrisk_stroke, cogdrisk_tbi, cogdrisk_dep, cogdrisk_stroke, 
         cogdrisk_afib, cogdrisk_smk), as_factor
  ) %>%
  select(NACCID, starts_with('cogdrisk'), z_cogdrisk)
  
## Combine datasets 
nacc_out <- nacc %>%
  left_join(nacc_lv, by = c('NACCID' = 'NACCID')) %>%
  left_join(nacc_lv2, by = c('NACCID' = 'NACCID')) %>%
  left_join(caide, by = 'NACCID') %>%
  left_join(mcaide, by = 'NACCID') %>%
  left_join(cogdrisk, by = 'NACCID')

## Export
write_csv(nacc_out, 'data/nacc.csv')
write_rds(nacc_out, 'data/nacc.rds', compress = 'gz')




































