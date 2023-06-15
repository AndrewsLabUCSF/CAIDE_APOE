library(tidyverse)
library(conflicted)
library(tabulizer)
library(glue)

farrer1997_path <- "~/Downloads/jama_278_16_041.pdf"
farrer1997_raw <- extract_tables(farrer1997_path, pages = 1)


test <- farrer1997_raw[[1]] %>%
  as_tibble() %>%
  slice(1:51) %>%
  unite(study, c(V1, V2, V3), sep = " ") %>%
  select(-V4, -V10) %>%
  unite(n_case, c(V5, V6), sep = "") %>%
  unite(race, c(V11:V15), sep = " ") %>%
  write_csv("~/Downloads/farrer1997_cohorts.csv")

farrer_raw <- read_csv('~/Dropbox/Research/UCSF/ABA/data/farrer1997_cohorts.csv')

farrer <- farrer_raw %>%
  mutate(
    n = cases + ctrls, 
    n_caucasian = n * caucasian, 
    n_aa = n * african_american,
    n_his = n * hispanic, 
    n_jap = n * japanese
    )


test <- farrer %>%
  group_by(ascertainment) %>%
  summarise(
    n = sum(n), 
    Caucasian = sum(n_caucasian), 
    'African American' = sum(n_aa), 
    Hispanic = sum(n_his), 
    Japanese = sum(n_jap), 
  ) %>%
  dplyr::filter(!is.na(ascertainment))

test %>% 
  gtsummary::tbl_summary(
    by = ascertainment
  ) 

tbl <- test %>%
  pivot_longer(., cols = c(n, Caucasian, 'African American', Hispanic, Japanese), 
             names_to = 'ancestry', values_to = 'n') %>%
  pivot_wider(names_from = ascertainment, values_from = n) %>%
  janitor::clean_names() %>%
  mutate(
    n = a + c + c_p + p, 
    a = round(a), 
    c = round(c), 
    p = round(p), 
    c_p = round(c_p),
    a_perc = round((a / n) * 100, 2),
    c_perc = round((c / n) * 100, 2),
    c_p_perc = round((c_p / n) * 100, 2),
    p_perc = round((p / n) * 100, 2), 
    Autopsy = glue("{a} ({a_perc}%)"),
    Clinical = glue("{c} ({c_perc}%)"),
    Population = glue("{p} ({p_perc}%)"),
    "Clinical/Population" = glue("{c_p} ({c_p_perc}%)")
  ) %>%
  select(ancestry, Autopsy, Clinical, Population, `Clinical/Population`)

gt::gt(tbl) %>%
  gt::gtsave("~/Downloads/ascertainment.png")























