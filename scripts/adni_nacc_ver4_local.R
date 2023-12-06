##########################################################
# Program to produce the combined ADNI and NACC dataset
##########################################################

library(tidyverse)
library(dplyr)
library(broom)
library(MASS)
library(gtsummary)
library(ggplot2)

setwd("/Users/elinorvelasquez/Desktop/adni_nacc/")
# setwd("/wynton/group/andrews/users/evelasquez/CAIDE_APOE/scripts/ADNI")

# cd ~/ for softlink
# (base) [evelasquez@log1 ~]$ ln -s /wynton/group/andrews/data/ resources
# adrc.raw <- read_csv("~/resources/NACC/data/ADRC_all_9_16_2022.csv") %>% 
#   janitor::clean_names()

nacc_data_raw <- read_csv("~/Desktop/adni_nacc/nacc.csv") 

adni_data_raw <- read_csv("~/Desktop/ADNI/adni_caide_calc.csv")

adni_novel_raw <- read_csv("~/Desktop/ADNI/adni_novel.csv")

adni_novel_raw$rid = as.factor(adni_novel_raw$rid)
adni_data_raw$rid = as.factor(adni_data_raw$rid)

adni_novel0 <- adni_novel_raw %>% dplyr::filter(viscode2 == "bl") 
                
adni_novel0 %>% count(diag) # adni_novel is new diagnosis
# 0   894
# 1  1514

adni_novel <- adni_data_raw %>%  left_join(adni_novel0, by=c('rid')) %>% 
  dplyr::filter(!is.na(diag))

adni_novel %>% count(diag)
# 0   873
# 1  1486

selected_adni <- adni_novel %>% 
  mutate(cohort = "ADNI") %>%
  dplyr::select( 
    rid, apoe, race, 
    i_age, i_education, gender,
    i_bmi, i_systbp, diag,
    hypchol3, hypchol2, cohort, 
    caide3, caide1, caide2,
    caide_apoe3, caide_apoe1, caide_apoe2,
    z_caide3, z_caide1, z_caide2, 
    mcaide3, mcaide1, mcaide2, 
    mcaide_apoe3, mcaide_apoe1, mcaide_apoe2,
    z_mcaide3, z_mcaide1, z_mcaide2 
  ) %>% 
  janitor::clean_names()

str(selected_adni) # 2369 x 30

selected_nacc <- nacc_data_raw %>% dplyr::filter(NACCVNUM == '1') %>%
	mutate(
	  dx_bl = case_when(
	    NACCUDSD == 1 ~ "CN",
	    NACCUDSD == 2 ~ "Impaired-not-MCI",
	    NACCUDSD == 3 ~ "MCI",
	    NACCUDSD == 4 ~ "ADRD",
	    TRUE ~ NA_character_
	  ),
	 diag = case_when(
	        dx_bl == 'AD' ~ 1,
	        dx_bl == 'ADRD' ~ 1,
	        dx_bl == 'CN' ~ 0,
	        dx_bl == 'EMCI' ~ 0,
	        dx_bl == 'LMCI' ~ 1,
	        dx_bl == "MCI" ~ 1,
	        dx_bl == 'SMC' ~ 0
	  ),
		cohort = "NACC",
		gender = case_when(
		  SEX == 1 ~ "Male",
		  SEX == 2 ~ "Female"
		),
		hypchol3 = i_hypchol,
		hypchol2 = i_hypchol,
		caide3 = caide,
		caide2 = caide,
		caide1 = caide,
		mcaide3 = mcaide,
		mcaide2 = mcaide,
		mcaide1 = mcaide,
		z_caide3 = z_caide,
		z_caide2 = z_caide,
		z_caide1 = z_caide,
		z_mcaide3 = z_mcaide,
		z_mcaide2 = z_mcaide,
		z_mcaide1 = z_mcaide,
		caideXapoe = case_when(
		  caide_apoe == "high_e2+" ~ "High, e2+",
		  caide_apoe == "high_e3/e3" ~ "High, e3/e3",
		  caide_apoe == "high_e4+" ~ "High, e4+",
		  caide_apoe == "low_e2+" ~ "Low, e2+",
		  caide_apoe == "low_e3/e3" ~ "Low, e3/e3",
		  caide_apoe == "low_e4+" ~ "Low, e4+",
		  caide_apoe == "mid_e2+" ~ "Mid, e2+",
		  caide_apoe == "mid_e3/e3" ~ "Mid, e3/e3",
		  caide_apoe == "mid_e4+" ~ "Mid, e4+",
		),
		mcaideXapoe = case_when(
		  mcaide_apoe == "high_e2+" ~ "High, e2+",
		  mcaide_apoe == "high_e3/e3" ~ "High, e3/e3",
		  mcaide_apoe == "high_e4+" ~ "High, e4+",
		  mcaide_apoe == "low_e2+" ~ "Low, e2+",
		  mcaide_apoe == "low_e3/e3" ~ "Low, e3/e3",
		  mcaide_apoe == "low_e4+" ~ "Low, e4+",
		  mcaide_apoe == "mid_e2+" ~ "Mid, e2+",
		  mcaide_apoe == "mid_e3/e3" ~ "Mid, e3/e3",
		  mcaide_apoe == "mid_e4+" ~ "Mid, e4+",
		),
		  caide_apoe3 = caideXapoe,
		  caide_apoe2 = caideXapoe,
		  caide_apoe1 = caideXapoe,
		  mcaide_apoe3 = mcaideXapoe,
		  mcaide_apoe2 = mcaideXapoe,
		  mcaide_apoe1 = mcaideXapoe,
  ) %>%
  dplyr::select(-c("caide_apoe")) %>%
  dplyr::select(-c("mcaide_apoe")) %>%
  dplyr::rename( 
      rid = NACCID, 
      i_age = NACCAGE, i_education = i_EDUC,
      i_bmi = i_NACCBMI, i_systbp = i_BPSYS
  ) %>%
  dplyr::select(-c(HYPCHOL)) %>%
  dplyr::select(
    rid, apoe, race, 
    i_age, i_education, gender,
    i_bmi, i_systbp, diag,
    hypchol3, hypchol2, cohort, 
    caide3, caide1, caide2,
    caide_apoe3, caide_apoe2, caide_apoe1,
    z_caide3, z_caide2, z_caide1, 
    mcaide3, mcaide1, mcaide2, 
    mcaide_apoe3, mcaide_apoe2, mcaide_apoe1,
    z_mcaide3, z_mcaide2, z_mcaide1 
  ) %>%
  janitor::clean_names()

str(selected_nacc) # 18237 x 30

adni_nacc <- rbind(selected_adni, selected_nacc)
write_csv(adni_nacc, 
          "~/Desktop/adni_nacc/adni_nacc.csv")

selected_adni %>% count(cohort) 
selected_nacc %>% count(cohort) 
adni_nacc %>% count(race) 




