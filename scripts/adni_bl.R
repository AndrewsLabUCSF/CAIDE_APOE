
library("tidyverse")
library("missForest")
library("glue")
library("dplyr")
library("broom")
library("ggplot2")
library("gtsummary")
library("rmarkdown")
library(lubridate)
library(vcfR)

setwd("~/gitcode/CAIDE_APOE")

## Import datasets
### Key variables 
merge.raw <- read_csv("resources/ADNI/ADNIMERGE_25Sep2023.csv") %>%
  janitor::clean_names()

### Diagnostic summary
dxsum_raw <- read_csv("resources/ADNI/DXSUM_PDXCONV_ADNIALL_25Sep2023.csv", guess_max = 14000) %>%
  janitor::clean_names() 

### patient demographics
ptdemo.raw <- read_csv("resources/ADNI/PTDEMOG_21Nov2023.csv") %>%
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
  select(phase, rid, ptid, viscode, apoe_geno) %>%
  mutate(apoe = case_when(
    apoe_geno %in% c(33) ~ 'e3/e3', 
    apoe_geno %in% c(24, 42, 43, 34, 44) ~ "e4+", 
    apoe_geno %in% c(23, 32, 22) ~ "e2+"
  ))

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

## CI Age of Onset
# ptdemo.raw %>% count(viscode) %>% arrange(-n) %>% print(n = Inf)
aoo <- ptdemo.raw %>% 
  group_by(ptid) %>%
  slice_min(visdate) %>%
  ungroup() %>%
  mutate(
    ptadbeg = na_if(ptadbeg, -4), 
    ptadbeg = na_if(ptadbeg, -1), 
    ptaddx = na_if(ptaddx, 9999), 
    ptdobyy = as.character(ptdobyy),
    ptdobyy = as.Date(ptdobyy, format = "%Y"),
    age = year(visdate) - year(ptdobyy),
    ad_yoo = case_when(
      !is.na(ptadbeg) ~ ptadbeg, 
      !is.na(ptaddx) ~ ptaddx
    ), 
    ad_yoo = as.character(ad_yoo),
    ad_yoo = as.Date(ad_yoo, format = "%Y"),
    ad_aoo = year(ad_yoo) - year(ptdobyy)
    ) %>%
  select(phase, ptid, viscode2, visdate, ptdobyy, age, ptcogbeg, ptmcibeg, ptadbeg, ptaddx, ad_yoo, ad_aoo
         ) 

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

### Vitals 
vitals <- vitals.raw %>%
  select(phase, rid, viscode2, visdate, vsweight, vswtunit, vsheight, vshtunit, vsbpsys, vsbpdia) %>%
  arrange(rid, visdate) %>%
  mutate(
    vsweight = na_if(vsweight, -1),
    vsweight = na_if(vsweight, -4),
    vsheight = na_if(vsheight, -1),
    vsheight = na_if(vsheight, -4),   
    vswtunit = na_if(vswtunit, -1),
    vshtunit = na_if(vshtunit, -1),
    vswtunit = na_if(vswtunit, -4),
    vshtunit = na_if(vshtunit, -4),
    vsbpsys = na_if(vsbpsys, -1),
  ) %>%
  group_by(rid) %>%
  tidyr::fill(vsheight, .direction = "down") %>%
  tidyr::fill(vshtunit, .direction = "down") %>%
  ungroup() %>%
  filter(viscode2 == "bl") %>%
  ## Fix stupid coding errors 
  mutate(
    vswtunit = case_when(
      rid == 6049 ~ 1,
      TRUE ~ vswtunit 
    ), 
    vshtunit = case_when(
      rid == 7123 ~ 1, 
      rid == 6916 ~ 2, 
      rid == 6134 ~ 2,
      TRUE ~ vshtunit 
    ),
  ) %>%
  mutate(
    weight_in_kg = case_when(vswtunit == 1 ~ vsweight*0.453592, #1 is lbs
                             vswtunit == 2 ~ vsweight,
                             TRUE ~ NA_real_
    ),
    height_in_cm = case_when(vshtunit == 1 ~ vsheight*2.54, #1 is inches
                             vshtunit == 2 ~ vsheight,
                             TRUE ~ NA_real_
    ),
    bmi = weight_in_kg / (height_in_cm / 100)^2
  )


### Cholesterol 
chol <- labs.raw %>%
  select(rid, phase, viscode, viscode2, examdate, rct20) %>%
  arrange(rid, examdate) %>%
  mutate(
    rct20 = as.numeric(rct20), 
    rct20 = na_if(rct20, -1),
    total_c = rct20 * 0.0259, 
    high_chol = ifelse(rct20 >= 240, 1, 0)
  ) %>%
  filter(!is.na(rct20)) %>%
  group_by(rid) %>%
  slice(which.min(examdate)) %>%
  ungroup()

## Merged data 
dxsum_merge <- dxsum %>% filter(viscode2 == "bl") %>% select(ptid, viscode2, naccudsd, naccetpr)
apoe_merge <- select(apoeres, -rid, -phase, -viscode) 
vitals_merge <- select(vitals, rid, vsbpsys, vsbpdia, weight_in_kg, height_in_cm, bmi)
chol_merge <- select(chol, rid, total_c, high_chol)
aoo_merge <- select(aoo, ptid, ad_aoo)

adni_wrangle <- merge %>% 
  left_join(dxsum_merge, by = "ptid") %>%
  left_join(aoo_merge, by = "ptid") %>%
  left_join(apoe_merge, by = "ptid") %>%
  left_join(vitals_merge, by = "rid") %>%
  left_join(chol_merge, by = "rid") %>%
  filter(race != "Other")  %>%
  filter(age >= 55 & !is.na(age)) %>%
  filter(is.na(ad_aoo) | ad_aoo >= 55) %>%
  filter(!is.na(apoe_geno)) %>%
  filter(naccetpr %in% c(1, 2, 4, 5, 6, 7, 8, 88))


### impute missing data 
require(doParallel)
registerDoParallel(cores=3) 
getDoParWorkers()
require(doRNG)
registerDoRNG(seed = 333)

adni_ximp <- adni_wrangle %>% 
  select(age, race, cdrsb, apoe_geno, ptgender, pteducat, height_in_cm, weight_in_kg, bmi, vsbpsys, vsbpdia, total_c) %>%
  mutate(
    race = as_factor(race), 
    apoe_geno = as_factor(apoe_geno), 
    ptgender = as_factor(ptgender), 
  ) %>%
  as.data.frame() %>%
  missForest(xmis = ., variablewise = TRUE, verbose = TRUE, xtrue = filter(., complete.cases(.)), parallelize = 'forests')

adni_imp <- adni_ximp %>%
  magrittr::extract2(1) %>%
  select(-age, -race, -cdrsb, -apoe_geno, -ptgender) %>%
  magrittr::set_colnames(paste0("i_", colnames(.))) %>%
  as_tibble() %>%
  mutate(
    i_pteducat = round(i_pteducat), 
    i_bmi = round(i_bmi, 1)
  )

adni <- adni_wrangle %>%
  bind_cols(adni_imp)

### CAIDE risk score
caide <- adni %>%
  mutate(
    caide_age = case_when(
      age <= 46 ~ 0, 
      age >= 47 & age < 53 ~ 3, 
      age >= 53 ~ 4, 
      TRUE ~ NA_real_
    ), 
    caide_educ = case_when(
      i_pteducat <= 6 ~ 3, 
      i_pteducat >= 7 & i_pteducat <= 9 ~ 2, 
      i_pteducat > 9 ~ 0, 
      TRUE ~ NA_real_
    ), 
    caide_sex = case_when(
      ptgender == "Male" ~ 1, 
      ptgender == "Female" ~ 0, 
      TRUE ~ NA_real_
    ), 
    caide_bmi = case_when(
      i_bmi < 30 ~ 0, 
      i_bmi >= 30 ~ 2, 
      TRUE ~ NA_real_
    ),
    caide_sbp = case_when(
      i_vsbpsys < 140 ~ 0, 
      i_vsbpsys >= 140 ~ 2, 
      TRUE ~ NA_real_
    ),
    i_caide_chol = case_when(
      i_total_c < 6.21 ~ 0, 
      i_total_c >= 6.21 ~ 2, 
      TRUE ~ NA_real_
    ), 
    caide_chol = case_when(
      total_c < 6.21 ~ 0, 
      total_c >= 6.21 ~ 2, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    caide = sum(caide_age, caide_educ, caide_sex, caide_bmi, caide_sbp, i_caide_chol, na.rm = F), 
    caide_u_chol = sum(caide_age, caide_educ, caide_sex, caide_bmi, caide_sbp, caide_chol, na.rm = F), 
    caide_wo_chol = sum(caide_age, caide_educ, caide_sex, caide_bmi, caide_sbp, na.rm = F), 
    caide_missing = sum(is.na(age), is.na(pteducat), is.na(ptgender), 
                        is.na(bmi), is.na(total_c), is.na(vsbpsys)),
  ) %>%
  ungroup() %>%
  select(ptid, i_caide_chol, starts_with("caide"))

## mCAIDE
mcaide <- adni %>%
  mutate(
    mcaide_age = case_when(
      age <= 64.5 ~ 0, 
      age > 64.5 & age < 73 ~ 1, 
      age >= 73 ~ 2, 
      TRUE ~ NA_real_
    ), 
    mcaide_educ = case_when(
      i_pteducat <= 11 ~ 2, 
      i_pteducat >= 12 & i_pteducat <= 16 ~ 1, 
      i_pteducat > 16 ~ 0, 
      TRUE ~ NA_real_
    ), 
    mcaide_sex = case_when(
      ptgender == "Male" ~ 1, 
      ptgender == "Female" ~ 0, 
      TRUE ~ NA_real_
    ), 
    mcaide_bmi = case_when(
      i_bmi < 30 ~ 0, 
      i_bmi >= 30 ~ 2, 
      TRUE ~ NA_real_
    ),
    mcaide_sbp = case_when(
      i_vsbpsys < 140 ~ 0, 
      i_vsbpsys >= 140 ~ 2, 
      TRUE ~ NA_real_
    ),
    i_mcaide_chol = case_when(
      i_total_c < 6.21 ~ 0, 
      i_total_c >= 6.21 ~ 2, 
      TRUE ~ NA_real_
    ), 
    mcaide_chol = case_when(
      total_c < 6.21 ~ 0, 
      total_c >= 6.21 ~ 2, 
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    mcaide = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_bmi, mcaide_sbp, i_mcaide_chol, na.rm = F), 
    mcaide_u_chol = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_bmi, mcaide_sbp, mcaide_chol, na.rm = F), 
    mcaide_wo_chol = sum(mcaide_age, mcaide_educ, mcaide_sex, mcaide_bmi, mcaide_sbp, na.rm = F), 
    mcaide_missing = sum(is.na(age), is.na(pteducat), is.na(ptgender), 
                        is.na(bmi), is.na(total_c), is.na(vsbpsys)),
  ) %>%
  ungroup() %>%
  select(ptid, i_mcaide_chol, starts_with("mcaide"))

## export datsets 
adni_out <- adni %>%
  left_join(caide, by = "ptid") %>%
  left_join(mcaide, by = "ptid") %>%
  mutate(
    htn = case_when(
      vsbpsys < 140 ~ FALSE, 
      vsbpsys >= 140 ~ TRUE, 
      TRUE ~ NA_real_
    ), 
    hld = case_when(
      total_c < 6.21 ~ FALSE, 
      total_c >= 6.21 ~ TRUE, 
      TRUE ~ NA_real_
    ), 
    i_htn = case_when(
      i_vsbpsys < 140 ~ FALSE, 
      i_vsbpsys >= 140 ~ TRUE, 
      TRUE ~ NA_real_
    ), 
    i_hld = case_when(
      i_total_c < 6.21 ~ FALSE, 
      i_total_c >= 6.21 ~ TRUE, 
      TRUE ~ NA_real_
    ), 
    dx = case_when(
      naccudsd == 1 ~ 0, 
      naccudsd == 3 ~ 1,
      naccudsd == 4 ~ 1
    )
  )
    

write_tsv(adni_out, "data/adni.csv")

## Table 1
adni_out %>%
  mutate(
    race = as_factor(race), 
    race = fct_recode(race, "NHW" = "Non-Hispanic White"), 
    apoe = as_factor(apoe), 
    apoe = fct_relevel(apoe, "e3/e3"), 
    obesity = case_when(
      bmi < 30 ~ FALSE, 
      bmi >= 30 ~ TRUE, 
      TRUE ~ NA_real_
    ), 
    dx = case_when(
      naccudsd == 1 ~ "CU", 
      naccudsd == 3 ~ "MCI",
      naccudsd == 4 ~ "ADRD"
    ), 
    dx = fct_relevel(dx, "CU", "MCI", "ADRD"),
    adrd = case_when(
      naccetpr == 88 ~ "CU", 
      naccetpr == 1 ~ "AD", 
      naccetpr == 4 ~ "PSP", 
      naccetpr == 7 ~ "FTLD", 
      naccetpr == 8 ~ "VCID", 
      TRUE ~ NA_character_
    ), 
    adrd = fct_relevel(adrd, "CU"), 
  ) %>%
  dplyr::select(race, age, ptgender, pteducat, htn, obesity, hld, caide, mcaide, apoe, dx, adrd) %>%
  tbl_summary(.,
              by = race,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"
              ), 
              label = list(ptgender ~ "Gender",
                           age ~ "Age",
                           pteducat ~ "Education",
                           race ~ "Race",
                           dx ~ "Diagnosis",
                           adrd ~ "ADRD",
                           obesity ~ "Obesity",
                           htn ~ "Hypertension",
                           hld ~ "Dyslipidemia",
                           caide ~ "CAIDE",
                           mcaide ~ "mCAIDE",
                           apoe ~ "APOE"
                           )
  )




