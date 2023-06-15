library(tidyverse)
library(janitor)

`%nin%` = negate(`%in%`)

setwd('~/Dropbox/Research/UCSF/ABA')

# NACC
## Import datasets
nacc.raw <- read_csv('resources/investigator_nacc58.csv')

## NACC variables of interest
nacc_df_lv <- select(nacc.raw, NACCID, NACCADC, VISITYR, PACKET, NACCVNUM,
               NACCAGE, 
               NACCMMSE, CDRSUM,
               NACCUDSD, # Prevelant
               DECAGE, NACCETPR
) 

dat_wrangle <- nacc_df_lv %>%
  mutate(
    ## Missing data 
    NACCMMSE = ifelse(between(NACCMMSE, 0, 30), NACCMMSE, NA), 
    DECAGE = ifelse(DECAGE %in% c(888, 999), NA, DECAGE),
  ) %>%
  ## Retain only information from last visit
  filter(NACCETPR %nin% c(99)) %>%
  group_by(NACCID) %>%
  filter(NACCVNUM == which.max(NACCVNUM)) %>%
  ungroup() %>%
  ## Age of Onset, for diagnosis, age at which symptoms began, for controls, age at last examination. 
  mutate(
    aoo = case_when(
      NACCETPR %in% c(1:30) ~ DECAGE, 
      NACCETPR == 88 ~ NACCAGE, 
      TRUE ~ NA_real_
    ), 
    dx = ifelse(NACCUDSD %in% c(3,4), 1, 0),
  ) 


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
exclude_ID <- filter(test, (NACCVNUM == 1 & NACCUDSD %in% c(3, 4)) | ADRD == 88) %>% distinct(NACCID) %>% pull(NACCID)
exclude_ADRD_ID <- filter(test, ADRD %in% c(1, 88)) %>% distinct(NACCID) %>% pull(NACCID)

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


ggplot(test, aes(x = aoo, y = as.factor(dx_aoo),  fill = as.factor(dx_aoo))) + 
  geom_violin()

ggplot() + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_smooth(data = nacc_lv %>% filter(dx_lv == 0), aes(x = NACCAGE_lv, y = ad_age_lv), color = 'steelblue') +
  geom_smooth(data = nacc_lv %>% filter(dx_lv == 1), aes(x = NACCAGE_lv, y = ad_age_lv), color = 'firebrick') +
  theme_bw()




































