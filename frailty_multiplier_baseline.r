##################################################################################################
############# Aim: derive the multiplier for the baseline frailty parameter

rm(list=ls())

library(tidyr)
library(dplyr)
library(HMDHFDplus)
library(stringr)
library(purrr)
library(zoo)
library(LifeTables)

data.dir <-"your directory"
out.dir <-"your directory"
plot.dir <-"your directory"

setwd(data.dir)

ourcohorts <- 1885:1954
years <- 1980:2019 
ages <- 65:95


# read the data from HMD
username <-  # your username
password <- # f your password



###########################cohort lt
## estimate cohort life tables from data for mortality rates for cohorts - no life tables in HMD for cohorts older than 1930
## and a0 from period from HMD
countries <- c("CHE", "DNK", "FIN", "FRATNP", "GBRTENW", "ITA","NLD","NOR","SWE")  


i=1
countryi <- countries[i]

lx_coh <- readHMDweb(CNTRY=countryi, item="cMx_1x1", username=username, password=password) %>% #read cohort mx data
  mutate(fem_M=as.numeric(Female), male_M=as.numeric(Male)) %>%
  select(Year,Age,fem_M,male_M) %>%
  mutate(Age= ifelse(Age=="110+","110",Age)) %>%
  filter(Year %in% ourcohorts) %>%
  pivot_longer(cols = c(fem_M, male_M),
               names_to = "sex",
               names_pattern = "(fem|male)_M",
               values_to = "mx") %>%
  #read in period a0 data
  left_join(readHMDweb(CNTRY=countryi, item="fltper_1x1", username=username, password=password) %>%
            filter(Age==0, Year %in% ourcohorts) %>%
            mutate(fem_a0=ax) %>%
            select(fem_a0, Year)%>%
  left_join(readHMDweb(CNTRY=countryi, item="mltper_1x1", username=username, password=password) %>%
              filter(Age==0, Year %in% ourcohorts) %>%
              mutate(male_a0=ax) %>%
              select(male_a0, Year)) %>%
     pivot_longer(cols = c(fem_a0, male_a0),
                 names_to = "sex",
                 names_pattern = "(fem|male)_a0",
                 values_to = "ax")) %>%
  mutate(ax=ifelse(Age==0,ax,0.5)) %>%
  mutate(mx = ifelse(is.na(mx), 0, mx)) %>%
  group_by(sex, Year) %>%
  group_split() %>%
  map_dfr(function(df) {
    out <- LifeTables::lt.mx(nmx = df$mx, nax = df$ax, age = df$Age)
    lt <- out$lt
    tibble(Year = unique(df$Year), sex  = unique(df$sex), Age=lt[, 1],  lx= lt[, 7])}) %>%
  mutate(country=countryi, cohort=Year, Year=Year+Age, lx_c=lx/100000) %>%
  filter(Year %in% years, Age %in%ages) %>%
  dplyr::select(country,sex,Year,Age,lx_c,cohort)

write.table(lx_coh, file="cohlt.csv",sep=",", row.names=FALSE)

for (i in 2:length(countries)){
  countryi <- countries[i]
  
  lx_coh <- readHMDweb(CNTRY=countryi, item="cMx_1x1", username=username, password=password) %>% #read cohort mx data
    mutate(fem_M=as.numeric(Female), male_M=as.numeric(Male)) %>%
    select(Year,Age,fem_M,male_M) %>%
    mutate(Age= ifelse(Age=="110+","110",Age)) %>%
    filter(Year %in% ourcohorts) %>%
    pivot_longer(cols = c(fem_M, male_M),
                 names_to = "sex",
                 names_pattern = "(fem|male)_M",
                 values_to = "mx") %>%
    #read in period a0 data
    left_join(readHMDweb(CNTRY=countryi, item="fltper_1x1", username=username, password=password) %>%
                filter(Age==0, Year %in% ourcohorts) %>%
                mutate(fem_a0=ax) %>%
                select(fem_a0, Year)%>%
                left_join(readHMDweb(CNTRY=countryi, item="mltper_1x1", username=username, password=password) %>%
                            filter(Age==0, Year %in% ourcohorts) %>%
                            mutate(male_a0=ax) %>%
                            select(male_a0, Year)) %>%
                pivot_longer(cols = c(fem_a0, male_a0),
                             names_to = "sex",
                             names_pattern = "(fem|male)_a0",
                             values_to = "ax")) %>%
    mutate(ax=ifelse(Age==0,ax,0.5)) %>%
    mutate(mx = ifelse(is.na(mx), 0, mx)) %>%
    group_by(sex, Year) %>%
    group_split() %>%
    map_dfr(function(df) {
      out <- LifeTables::lt.mx(nmx = df$mx, nax = df$ax, age = df$Age)
      lt <- out$lt
      tibble(Year = unique(df$Year), sex  = unique(df$sex), Age=lt[, 1],  lx= lt[, 7])}) %>%
    mutate(country=countryi, cohort=Year, Year=Year+Age, lx_c=lx/100000) %>%
    filter(Year %in% years, Age %in%ages) %>%
    dplyr::select(country,sex,Year,Age,lx_c,cohort)
  
  write.table(lx_coh, file="cohlt.csv",sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

#######################################################################################
############### multiplier for the baseline frailty

multip <- read.table(file="cohlt.csv",sep=",", header=TRUE) %>%
  mutate(Age=as.numeric(Age)) %>%
  group_by(country, sex, Year) %>%
  arrange(Age, .by_group = TRUE) %>%
  mutate(multip = (lx_c + lead(lx_c))/2) %>%   # last row in each group gets NA
  ungroup() %>%
  filter(Age<95) %>%
  select(-cohort)

setwd(data.dir)  
write.table(multip, file="multip.csv",sep=",",row.names=FALSE)  
  
  
  
  