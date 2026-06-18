###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Download HMD deaths and exposures by Lexis triangle for the nine study populations, restrict to the study window, reshape to a long
##           (country x sex x triangle) table, place each triangle at its mid-point age / period / cohort coordinate, and write out in one CSV.
## Inputs:   HMD items: "Exposures_lexis", "Deaths_lexis".
## Output:   <out.dir>/lexis_1960.csv. 
## Window:   Ages 65-94; periods 1960-2019;
## Depends:  HMDHFDplus, dplyr, tidyr
###############################################################################
rm(list = ls())

library(HMDHFDplus) 
library(dplyr) 
library(tidyr)

## ---- User settings ---------------------------------------------------------
out.dir  <- #"your/output/directory"   
username <-"your_HMD_username"        
password <- "your_HMD_password"

## ---- Study design ----------------------------------------------------------
## HMD codes -> CHE=Switzerland, DNK=Denmark, FIN=Finland, FRATNP=France (total pop.), GBRTENW=England & Wales, ITA=Italy, 
## NLD=Netherlands, NOR=Norway, SWE=Sweden.
countries  <- c("CHE","DNK","FIN","FRATNP","GBRTENW","ITA","NLD","NOR","SWE")
ages       <- 65:94       # 65 = old-age focus; 94 = data-quality upper bound
years      <- 1960:2019
## Cohort range follows from age x period: age 94 in 1960 -> cohort 1885; age 65 in 2019 -> cohort 1954. 
ourcohorts <- 1865:1954

## ---- Download + shape one country ------------------------------------------
## Each (Year, Age) cell holds two Lexis triangles, identified by Cohort:
##   Cohort == Year - Age      -> LOWER triangle (born in calendar year Year-Age)
##   Cohort == Year - Age - 1  -> UPPER triangle (born one year earlier)
##
## Triangle centroids (Carstensen 2007), so the continuous dAPC splines see the
## triangle rather than the integer cell:
##   LOWER: Age + 1/3, Period + 2/3   (=> Cohort centroid = Period - Age = (Y-A) + 1/3)
##   UPPER: Age + 2/3, Period + 1/3   (=> Cohort centroid = Period - Age = (Y-A) - 1/3)
## Cohort is DERIVED as Period - Age after shifting
get_country <- function(cc) {
  exp <- readHMDweb(CNTRY = cc, item = "Exposures_lexis",
                    username = username, password = password) %>%
    filter(Cohort %in% ourcohorts, Year %in% years, Age %in% ages) %>%
    dplyr::select(-Total, -OpenInterval) %>%
    pivot_longer(c(Female, Male), names_to = "sex", values_to = "Exp")
  
  dth <- readHMDweb(CNTRY = cc, item = "Deaths_lexis",
                    username = username, password = password) %>%
    filter(Cohort %in% ourcohorts, Year %in% years, Age %in% ages) %>%
    dplyr::select(-Total, -OpenInterval) %>%
    pivot_longer(c(Female, Male), names_to = "sex", values_to = "Dx")
  
  ## Join on the INTEGER cell keys, then shift to triangle centroids.
  exp %>%
    left_join(dth, by = c("Year", "Age", "Cohort", "sex")) %>%
    mutate(
      country = cc,
      older   = ifelse(Cohort == Year - Age, 1L, 0L),  # 1 = lower, 0 = upper
      Age     = ifelse(older == 1L, Age  + 1/3, Age  + 2/3),
      Year    = ifelse(older == 1L, Year + 2/3, Year + 1/3),
      Cohort  = Year - Age
    )
}

## ---- Run over all selected countries---------------------------------
lexis <- bind_rows(lapply(countries, get_country))

write.table(lexis, file = file.path(out.dir, "lexis_1960.csv"),
            sep = ",", row.names = FALSE)