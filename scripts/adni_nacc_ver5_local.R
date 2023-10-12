#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(broom)
library(MASS)
library(gtsummary)
library(ggplot2)

# setwd("/Users/elinorvelasquez/Desktop/adni_nacc/")
# setwd("/wynton/group/andrews/users/evelasquez/CAIDE_APOE/scripts/ADNI")

# cd ~/ for softlink
# (base) [evelasquez@log1 ~]$ ln -s /wynton/group/andrews/data/ resources
# adrc.raw <- read_csv("~/resources/NACC/data/ADRC_all_9_16_2022.csv") %>% 
#   janitor::clean_names()

# m/caide computed in all these datasets
nacc_data_raw <- read_csv("~/nacc.csv") # output from nacc_bl.R
mod_nacc_data_raw <- read_csv("~/no_chol_nacc.csv")  # modified nacc data; output from modified_nacc_bl.R

# select the caide-no-chol variables (caide, mcaide, caide_apoe, mcaide_apoe, z_caide, z_mcaide)
# from mod_nacc_data and recode them as caide1, mcaide1, etc; join these variables to the nacc_data
mod_nacc_data <- mod_nacc_data_raw %>% 
			dplyr::select(NACCID, caide, mcaide, caide_apoe, mcaide_apoe, z_caide, z_mcaide) %>%
			dplyr::rename(caide1 = caide, mcaide1 = mcaide, caide_apoe1 = caide_apoe,
				mcaide_apoe1 = mcaide_apoe, z_caide1 = z_caide, z_mcaide1 = z_mcaide) 

nacc_data <- nacc_data_raw %>% left_join(mod_nacc_data, by=c('NACCID'))

# nacc_data %>% dplyr::select(starts_with("i_")) %>% summary()
# nacc_data %>% count(is.na(i_EDUC))

adni_data <- read_csv("~/adni_caide_calc.csv") # output from adni_ver3.R

selected_adni <- adni_data %>% 
  mutate(cohort = "ADNI") %>%
  dplyr::select( 
    rid, apoe, race, 
    i_age, i_education, gender,
    i_bmi, i_systbp, dx_bl,
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

selected_nacc <- nacc_data %>% dplyr::filter(NACCVNUM == '1') %>%
	mutate(
	  dx_bl = case_when(
	    NACCUDSD == 1 ~ "CN",
	    NACCUDSD == 2 ~ "Impaired-not-MCI",
	    NACCUDSD == 3 ~ "MCI",
	    NACCUDSD == 4 ~ "ADRD",
	    TRUE ~ NA_character_
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
		mcaide3 = mcaide,
		mcaide2 = mcaide,
		z_caide3 = z_caide,
		z_caide2 = z_caide,
		z_mcaide3 = z_mcaide,
		z_mcaide2 = z_mcaide,
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
		mcaide_apoe3 = mcaideXapoe,
		mcaide_apoe2 = mcaideXapoe,
		caideXapoe1 = case_when(
                  caide_apoe1 == "high_e2+" ~ "High, e2+",
                  caide_apoe1 == "high_e3/e3" ~ "High, e3/e3",
                  caide_apoe1 == "high_e4+" ~ "High, e4+",
                  caide_apoe1 == "low_e2+" ~ "Low, e2+",
                  caide_apoe1 == "low_e3/e3" ~ "Low, e3/e3",
                  caide_apoe1 == "low_e4+" ~ "Low, e4+",
                  caide_apoe1 == "mid_e2+" ~ "Mid, e2+",
                  caide_apoe1 == "mid_e3/e3" ~ "Mid, e3/e3",
                  caide_apoe1 == "mid_e4+" ~ "Mid, e4+",
                ),
                mcaideXapoe1 = case_when(
                  mcaide_apoe1 == "high_e2+" ~ "High, e2+",
                  mcaide_apoe1 == "high_e3/e3" ~ "High, e3/e3",
                  mcaide_apoe1 == "high_e4+" ~ "High, e4+",
                  mcaide_apoe1 == "low_e2+" ~ "Low, e2+",
                  mcaide_apoe1 == "low_e3/e3" ~ "Low, e3/e3",
                  mcaide_apoe1 == "low_e4+" ~ "Low, e4+",
                  mcaide_apoe1 == "mid_e2+" ~ "Mid, e2+",
                  mcaide_apoe1 == "mid_e3/e3" ~ "Mid, e3/e3",
                  mcaide_apoe1 == "mid_e4+" ~ "Mid, e4+",
                )
  ) %>%
  dplyr::select(-c("caide_apoe", "caide_apoe1")) %>%
  dplyr::select(-c("mcaide_apoe", "mcaide_apoe1")) %>%
  dplyr::rename( 
        rid = NACCID, 
        i_age = NACCAGE, i_education = i_EDUC,
        i_bmi = i_NACCBMI, i_systbp = i_BPSYS,
	caide_apoe1 = caideXapoe1,
	mcaide_apoe1 = mcaideXapoe1	
  ) %>%
  dplyr::select(-c(HYPCHOL)) %>%
  dplyr::select(
    rid, apoe, race, 
    i_age, i_education, gender,
    i_bmi, i_systbp, dx_bl,
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
write_csv(adni_nacc, "~/adni_nacc.csv")




