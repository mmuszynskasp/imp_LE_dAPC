###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Fit the MAIN detrended Age-Period-Cohort model for every country x sex: B-splines with 10 knots per dimension and a
##           drift extracted with weights proportional to exposures. Save the age, period and cohort effects (with their model-based CIs), the
##           drift (with CI), and the nested-model deviance (ANOVA) table.
## Input:    <out.dir>/lexis_1960.csv  (from 1_prepare_data.R) 
## Outputs:  written to <models.dir>: ageeffectbs1960.csv, pereffectbs1960.csv, cohorteffectbs1960.csv, driftbs1960.csv, myanovabs1960.csv
## Depends:  Epi, dplyr 
###############################################################################
rm(list = ls())

library(Epi) 
library(dplyr)

data.dir    <- #lexis file is in this directory
models.dir  <- #<- "your/models/directory"   
## ---- Load prepared Lexis data ----------------------------------------------
lexis     <- read.table(file.path(data.dir, "lexis_1960.csv"), sep = ",", header = TRUE)
countries <- unique(lexis$country)
sexes     <- unique(lexis$sex)

## ---- Fit one model per country x sex, accumulate, write once ---------------
age_list <- per_list <- coh_list <- drift_list <- anova_list <- fits <- list() 
k <- 0L

for (cc in countries) {
  for (sx in sexes) {
    k <- k + 1L
    
    ## apc.fit wants a plain data.frame with columns A, P, D, Y.
    ## D = deaths, Y = exposures (offset). Rounded to integers for the Poisson fit. 
    dat <- lexis %>%
      filter(country == cc, sex == sx) %>%
      transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>%
      as.data.frame()
    
    m <- apc.fit(data = dat,
                 npar = c(A = 10, P = 10, C = 10),
                 model = "bs", parm = "APC",
                 dr.extr = "Y", ref.p = 1960 + 1/3)
    
    fits[[paste(cc, sx, sep = "_")]] <- m     #store full fit data
    
    ## Keep ALL columns of each effect matrix: tabulation point, point estimate, and the 2.5% / 97.5% model-based CIs. 
    age_list[[k]]   <- data.frame(m$Age, country = cc, sex = sx, check.names = FALSE)
    per_list[[k]]   <- data.frame(m$Per, country = cc, sex = sx, check.names = FALSE)
    coh_list[[k]]   <- data.frame(m$Coh, country = cc, sex = sx, check.names = FALSE)
    
    ## Drift: keep estimate AND its CI (was Drift[1,1] only). 
    drift_list[[k]] <- data.frame(as.list(m$Drift[1, ]),
                                  country = cc, sex = sx, check.names = FALSE)
    
    ## Full nested-model ANOVA: deviances AND the LRT df / p-values that back the "non-linear period and cohort effects significant at a<=0.01" statements
    ## and the deviance-based contribution decomposition (Table A1 / Figure A4).
    anova_list[[k]] <- data.frame(m$Anova, country = cc, sex = sx, check.names = FALSE)
  }
}

## ---- Write outputs ------------------------------------
## Fitted objects for the CI step: carries m$Model (coef/vcov), m$Knots, m$Ref.
saveRDS(fits, file.path(models.dir, "fits_bs1960.rds"))

write.table(bind_rows(age_list),   file.path(models.dir, "ageeffectbs1960.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(per_list),   file.path(models.dir, "pereffectbs1960.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(coh_list),   file.path(models.dir, "cohorteffectbs1960.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(drift_list), file.path(models.dir, "driftbs1960.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(anova_list), file.path(models.dir, "myanovabs1960.csv"),
            sep = ",", row.names = FALSE)