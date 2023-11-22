library(tidyverse)
library(janitor)

## Import datasets 
nacc.raw <- read_csv('data/nacc.csv')
adni.raw <- read_tsv('data/adni.csv')


## ADNI 
adni.raw %>% count(adrd)
adni.raw 
