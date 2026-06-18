###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Fit the dAPC model to BASELINE (frailty-adjusted) mortality for every country x sex. 
##           Baseline deaths are observed deaths rescaled to the standard individual (frailty mean 1 at age 0; variance s2, here 0.25); baseline exponent is −s2.
##           Model as for Total rates: B-spline 10-knot and weighted-drift specification as the total model in 02.
## Inputs:   <out.dir>/lexis.csv     (from 1_prepare_data.R), <models.dir>/multip.csv (from 3_frailty_multiplier.R)
## Outputs:  written to <models.dir>: ageeffectindiv10w.csv, pereffectindiv10w.csv, cohorteffectindiv10w.csv, driftindiv10w.csv, myanovaindiv10w.csv,
##           fits_indiv10w.rds  (fitted objects for the CI step: coef/vcov, knots)
## Depends:  Epi, dplyr
###############################################################################
rm(list = ls())

library(Epi)
library(dplyr)

out.dir    <-  #"your/output/directory"   # location of lexis.csv
models.dir <-  #"your/models/directory"   # multip.csv is in this directory; outputs go here

s2 <- 0.25   # assumed frailty variance (main text); baseline exponent is -s2

## ---- Build baseline (frailty-adjusted) triangle data -----------------------
## Join the cell-level multiplier (integer year x age) onto the Lexis triangles by floored coordinates: both triangles in a cell receive the same multiplier.
## A = Age, P = Year keep the continuous centroid coordinates, exactly as in the total model (02); cohort is derived inside apc.fit as P - A.
pc_data_si <- read.table(file.path(out.dir, "lexis.csv"),
                         sep = ",", header = TRUE) %>%
  mutate(year = floor(Year), age = floor(Age)) %>%
  left_join(
    read.table(file.path(models.dir, "multip.csv"), sep = ",", header = TRUE) %>%
      rename(year = Year, age = Age) %>%
      mutate(sex = ifelse(sex == "fem", "Female", "Male")),
    by = c("country", "sex", "year", "age")
  ) %>%
  mutate(Dxi = Dx * (multip^(-s2)),   # baseline deaths
         D   = round(Dxi),
         Y   = round(Exp),
         A   = Age,
         P   = Year)

## Guard: an unmatched cell gives multip = NA -> D = NA, which apc.fit would reject with an opaque error. Fail early with a clear message instead.
if (anyNA(pc_data_si$D))
  stop("NA baseline deaths after the multiplier join -- check year/age/sex matching.")

## ---- Fit one model per country x sex, accumulate ---------------
countries <- unique(pc_data_si$country)
sexes     <- unique(pc_data_si$sex)

age_list <- per_list <- coh_list <- drift_list <- anova_list <- fits <- list()
k <- 0L

for (cc in countries) {
  for (sx in sexes) {
    k <- k + 1L
    
    dat <- pc_data_si %>%
      filter(country == cc, sex == sx) %>%
      transmute(A, P, D, Y) %>%
      as.data.frame()
    
    m <- apc.fit(data = dat,
                 npar = c(A = 10, P = 10, C = 10),
                 model = "bs", parm = "APC",
                 dr.extr = "Y", ref.p = 1980 + 1/3)
    
    ## Keep all effect columns (tabulation point, estimate, 2.5% / 97.5% CIs).
    age_list[[k]]   <- data.frame(m$Age, country = cc, sex = sx, check.names = FALSE)
    per_list[[k]]   <- data.frame(m$Per, country = cc, sex = sx, check.names = FALSE)
    coh_list[[k]]   <- data.frame(m$Coh, country = cc, sex = sx, check.names = FALSE)
    drift_list[[k]] <- data.frame(as.list(m$Drift[1, ]),
                                  country = cc, sex = sx, check.names = FALSE)
    anova_list[[k]] <- data.frame(m$Anova, country = cc, sex = sx, check.names = FALSE)
    
    ## Fitted object for the CI step (coef/vcov via m$Model, plus m$Knots, m$Ref).
    fits[[paste(cc, sx, sep = "_")]] <- m
  }
}

## ---- Write outputs ---------------------------------------------------------
saveRDS(fits, file.path(models.dir, "fits_indiv10w.rds"))

write.table(bind_rows(age_list),   file.path(models.dir, "ageeffectindiv10w.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(per_list),   file.path(models.dir, "pereffectindiv10w.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(coh_list),   file.path(models.dir, "cohorteffectindiv10w.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(drift_list), file.path(models.dir, "driftindiv10w.csv"),
            sep = ",", row.names = FALSE)
write.table(bind_rows(anova_list), file.path(models.dir, "myanovaindiv10w.csv"),
            sep = ",", row.names = FALSE)
