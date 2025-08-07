rm(list=ls())

library(tidyr)
library(dplyr)
library(HMDHFDplus)

out.dir <- #your path
plot.dir <- #your path

username <- # replace with your username
password <- # replace with your password


countries <- c("CHE", "DNK", "FIN","FRATNP", "GBRTENW","ITA","NLD","NOR","SWE")  
sexes <- c("Female","Male")
  
ourcohorts <- 1885:1973
years <- 1980:2019 
ages <- 65:94


i=1
countrysel <- countries[i]
lexis <- readHMDweb(CNTRY = countrysel, item = "Exposures_lexis", username = username, password = password) %>%
  filter(Cohort %in% ourcohorts, Year %in% years, Age %in%ages) %>%
  dplyr::select(-c(Total, OpenInterval)) %>%
  pivot_longer(cols=Female:Male, values_to = "Exp", names_to = "sex") %>%
  mutate(country=countrysel) %>%
  left_join(readHMDweb(CNTRY = countrysel, item = "Deaths_lexis", username = username, password = password) %>%
               filter(Cohort %in% ourcohorts, Year %in% years, Age %in%ages) %>%
               dplyr::select(-c(Total, OpenInterval)) %>%
               pivot_longer(cols=Female:Male, values_to = "Dx", names_to = "sex") %>%
               mutate(country=countrysel)) %>%
  mutate(ageC=Year-Age,
         older=ifelse(Cohort==ageC,1,0),
         Year=ifelse(older==0,as.numeric(Year)+1/3,as.numeric(Year)+2/3),
         Age=ifelse(older==0,as.numeric(Age)+1/3,as.numeric(Age)+2/3),
         Cohort=ifelse(older==0,as.numeric(Cohort)+1/3,as.numeric(Cohort)+2/3))


setwd(out.dir)
write.table(lexis, file="lexis.csv", sep=",", row.names = FALSE)

for (i in 2:length(countries)){
countrysel <- countries[i]

lexis <- readHMDweb(CNTRY = countrysel, item = "Exposures_lexis", username = username, password = password) %>%
  filter(Cohort %in% ourcohorts, Year %in% years, Age %in%ages) %>%
  dplyr::select(-c(Total, OpenInterval)) %>%
  pivot_longer(cols=Female:Male, values_to = "Exp", names_to = "sex") %>%
  mutate(country=countrysel) %>%
  left_join(readHMDweb(CNTRY = countrysel, item = "Deaths_lexis", username = username, password = password) %>%
              filter(Cohort %in% ourcohorts, Year %in% years, Age %in%ages) %>%
              dplyr::select(-c(Total, OpenInterval)) %>%
              pivot_longer(cols=Female:Male, values_to = "Dx", names_to = "sex") %>%
              mutate(country=countrysel)) %>%
  mutate(ageC=Year-Age,
         older=ifelse(Cohort==ageC,1,0),
         Year=ifelse(older==0,as.numeric(Year)+1/3,as.numeric(Year)+2/3),
         Age=ifelse(older==0,as.numeric(Age)+1/3,as.numeric(Age)+2/3),
         Cohort=ifelse(older==0,as.numeric(Cohort)+1/3,as.numeric(Cohort)+2/3))
  write.table(lexis, file="lexis.csv", sep=",", row.names = FALSE, col.names=FALSE, append=TRUE)
}
