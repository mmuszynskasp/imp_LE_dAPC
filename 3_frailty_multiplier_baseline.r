###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Derive the survival term used to convert observed (total) mortality into baseline mortality under the gamma-frailty assumption. 
## Inputs:   HMD data. Items: cMx_1x1 (cohort rates), fltper_1x1 / mltper_1x1 (period life tables, used only for the age-0 ax). 
##           Cohort survival l(a,c) is built from HMD cohort rates, because HMD cohort life tables are unavailable for cohorts born before ~1930.
## Outputs:  written to <data.dir>: cohlt.csv  - cohort survival lx_c by country, sex, year, age, cohort
##           multip.csv - the averaged-survival multiplier by country, sex, year, age
## Window:   cohorts 1885-1954; years 1980-2019; ages 0-94.
## Depends:  HMDHFDplus, dplyr, tidyr, purrr, tibble, LifeTables
###############################################################################
rm(list = ls())

library(HMDHFDplus)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(LifeTables)

data.dir <- #"your/data/directory"     
username <- "your_HMD_username"        
password <- "your_HMD_password"


## ---- Study design ----------------------------------------------------------
ourcohorts <- 1885:1954
years      <- 1980:2019
ages       <- 65:95 # for the outcome survival, 0-94 for the multiplier
countries  <- c("CHE","DNK","FIN","FRATNP","GBRTENW","ITA","NLD","NOR","SWE")

## ---- Build one country's cohort survival -----------------------------------
## Steps:
##  1. cMx_1x1: cohort death rates, indexed by birth-cohort Year and Age.
##  2. Age-0 ax taken from the PERIOD life tables in the matching year; all other ages use ax = 0.5. (a0 affects only the level of lx accumulated from age 0.)
##  3. lt.mx builds a life table per (sex, cohort); lx is column 7 of its output.
##  4. cohort = birth year; calendar Year = cohort + Age; lx_c = lx / 100000 = S(a).

get_cohort_lt <- function(cc) {
  keep_years <- c(min(years) - 1L, years)   # 1979:2019, The extra prior year (1979) exists only so the same-age multiplier below is defined at the 1980 reference.
  
  mx <- readHMDweb(CNTRY = cc, item = "cMx_1x1",
                   username = username, password = password) %>%
    mutate(fem_M = as.numeric(Female), male_M = as.numeric(Male)) %>%
    select(Year, Age, fem_M, male_M) %>%
    mutate(Age = ifelse(Age == "110+", "110", Age)) %>%   
    filter(Year %in% ourcohorts) %>%
    pivot_longer(c(fem_M, male_M), names_to = "sex",
                 names_pattern = "(fem|male)_M", values_to = "mx")
  
  a0 <- readHMDweb(CNTRY = cc, item = "fltper_1x1",
                   username = username, password = password) %>%
    filter(Age == 0, Year %in% ourcohorts) %>%
    mutate(fem_a0 = ax) %>% select(fem_a0, Year) %>%
    left_join(readHMDweb(CNTRY = cc, item = "mltper_1x1",
                         username = username, password = password) %>%
                filter(Age == 0, Year %in% ourcohorts) %>%
                mutate(male_a0 = ax) %>% select(male_a0, Year),
              by = "Year") %>%
    pivot_longer(c(fem_a0, male_a0), names_to = "sex",
                 names_pattern = "(fem|male)_a0", values_to = "ax")
  
  mx %>%
    left_join(a0, by = c("Year", "sex")) %>%
    mutate(ax = ifelse(Age == 0, ax, 0.5),
           mx = ifelse(is.na(mx), 0, mx)) %>%
    group_by(sex, Year) %>%
    group_split() %>%
    map_dfr(function(df) {
      out <- LifeTables::lt.mx(nmx = df$mx, nax = df$ax, age = df$Age)
      lt  <- out$lt
      tibble(Year = unique(df$Year), sex = unique(df$sex),
             Age = lt[, 1], lx = lt[, 7])
    }) %>%
    mutate(country = cc, cohort = Year, Year = Year + Age, lx_c = lx / 100000) %>%
    filter(Year %in% keep_years, Age %in% ages) %>%
    dplyr::select(country, sex, Year, Age, lx_c, cohort)
}

##---- Run all countries -----------------------------------------------------
lx_coh <- map_dfr(countries, get_cohort_lt)

write.table(lx_coh %>% filter(Year %in% years),
            file.path(data.dir, "cohlt.csv"),
            sep = ",", row.names = FALSE)

## ---- Frailty multiplier ----------------------------------------------------
## Base of eq. (2): (l(a,c) + l(a,c+1))/2, with BOTH cohorts evaluated AT AGE a.
multip <- lx_coh %>%
  arrange(country, sex, Age, Year) %>%
  group_by(country, sex, Age) %>%
  mutate(multip = (lx_c + lag(lx_c)) / 2) %>%   # earliest year per group -> NA
  ungroup() %>%
  filter(Year %in% years, Age < 95) %>%
  dplyr::select(country, sex, Year, Age, lx_c, multip)

write.table(multip, file.path(data.dir, "multip.csv"),
            sep = ",", row.names = FALSE)
