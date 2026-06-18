###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Appendix table documenting the OBSERVED average annual change in
##           partial LE (ages 65-94) over each sub-period, paired with the
##           dAPC (model) Total for comparison -- backs the descriptive
##           acceleration/deceleration sentences in Section 3.1.
## Inputs:   <data.dir>/lexis.csv
##           <data.dir>/LE_change_decomp.csv   (for the dAPC Total column)
## Outputs:  <data.dir>/LE_observed_change.csv
##           <plot.dir>/obs_vs_dapc.tex   (single combined table; female + male columns)
## Depends:  dplyr, tidyr, kableExtra
###############################################################################
library(dplyr)
library(tidyr)
library(kableExtra)

  data.dir   <- #"your/data/directory"  
  plot.dir   <- #"your/plots/directory" 
  
cntry_map <- c(CHE = "Switzerland", DNK = "Denmark", FIN = "Finland",
               FRATNP = "France", GBRTENW = "England & Wales", ITA = "Italy",
               NLD = "Netherlands", NOR = "Norway", SWE = "Sweden")

periods <- data.frame(
  period = c("1983-1993", "1993-2003", "2003-2013", "2013-2017"),
  t0     = c(1983, 1993, 2003, 2013),
  t1     = c(1993, 2003, 2013, 2017),
  stringsAsFactors = FALSE
)
## sub-periods first, full-window summary last (matches the decomposition tables)
period_lv <- c(periods$period, "1983-2017")

lexis <- read.table(file.path(data.dir, "lexis.csv"), sep = ",", header = TRUE)

## ---- life-table helpers ----------------------------------------------------
## Partial LE over ascending ages, ax = 0.5, open last interval.
## (scalar form; verified identical to the matrix pLE_matrix used in script 05)
pLE <- function(m) {
  K  <- length(m)
  qx <- m / (1 + 0.5 * m); qx[K] <- 1
  px <- 1 - qx
  lx <- numeric(K); lx[1] <- 1
  if (K > 1) for (k in 2:K) lx[k] <- lx[k - 1] * px[k - 1]
  dx <- numeric(K)
  if (K > 1) dx[1:(K - 1)] <- lx[1:(K - 1)] - lx[2:K]
  dx[K] <- lx[K]
  Lx <- numeric(K)
  if (K > 1) Lx[1:(K - 1)] <- lx[2:K] + 0.5 * dx[1:(K - 1)]
  Lx[K] <- lx[K] / m[K]
  sum(Lx)
}

## observed partial LE by calendar year for one country x sex
obs_LE_by_year <- function(d) {
  ageint  <- floor(d$Age); yearint <- floor(d$Year)
  yrs <- sort(unique(yearint)); ags <- sort(unique(ageint))
  cellf <- factor(paste(yearint, ageint, sep = "_"))
  D <- rowsum(round(d$Dx),  cellf)[, 1]      # deaths per (year, age) cell
  Y <- rowsum(round(d$Exp), cellf)[, 1]      # exposure per cell
  rate <- D / Y
  LE <- vapply(yrs, function(y) pLE(rate[paste(y, ags, sep = "_")]), numeric(1))
  data.frame(year = yrs, LE = LE)
}

## average annual change in LE over [t0, t1]
pchg <- function(LE, years, t0, t1)
  (LE[match(t1, years)] - LE[match(t0, years)]) / (t1 - t0)

## ---- observed sub-period change for every country x sex --------------------
obs <- lexis %>%
  group_by(country, sex) %>%
  group_modify(~{
    le <- obs_LE_by_year(.x)
    sub <- data.frame(
      period   = periods$period,
      observed = mapply(function(a, b) pchg(le$LE, le$year, a, b),
                        periods$t0, periods$t1)
    )
    ## full-window 1983-2017 average annual change (telescopes to endpoint form)
    summ <- data.frame(period = "1983-2017",
                       observed = pchg(le$LE, le$year, 1983, 2017))
    bind_rows(sub, summ)
  }) %>%
  ungroup()

write.csv(obs, file.path(data.dir, "LE_observed_change.csv"), row.names = FALSE)

## ---- pair with the dAPC (model) Total --------------------------------------
dapc <- read.csv(file.path(data.dir, "LE_change_decomp.csv")) %>%
  filter(component == "total") %>%
  dplyr::select(country, sex, period, dapc = estimate)

tab <- obs %>%
  left_join(dapc, by = c("country", "sex", "period")) %>%
  mutate(country = recode(country, !!!as.list(cntry_map)),
         period  = factor(period, levels = period_lv),
         Observed       = sprintf("%.2f", observed),
         `dAPC (model)` = sprintf("%.2f", dapc)) %>%
  arrange(country, sex, period)

## ---- one combined table: female columns then male columns ------------------
## reshape so each (Country, Period) has Observed/dAPC for both sexes side by side
w <- tab %>%
  mutate(Period = gsub("-", "\u2013", as.character(period))) %>%
  dplyr::select(country, Period, period, sex, Observed, `dAPC (model)`) %>%
  pivot_wider(names_from = sex,
              values_from = c(Observed, `dAPC (model)`),
              names_glue = "{sex}_{.value}") %>%
  arrange(country, period) %>%
  dplyr::select(Country = country, Period,
                `Female_Observed`, `Female_dAPC (model)`,
                `Male_Observed`,   `Male_dAPC (model)`)

## column headers: drop the sex prefix; the sex grouping is added via add_header_above
colnames(w) <- c("Country", "Period",
                 "Observed", "dAPC (model)", "Observed", "dAPC (model)")

n_per_block <- length(period_lv)
summ_rows   <- seq(n_per_block, nrow(w), by = n_per_block)

cap <- paste(
  "Observed versus model (dAPC) average annual change in partial life",
  "expectancy (ages 65--94), by sex, in years per year.")

kbl(w, format = "latex", booktabs = TRUE, escape = TRUE, linesep = "",
    align = c("l", "l", "r", "r", "r", "r"),
    caption = cap,
    label = "obs_dapc") %>%
  kable_styling(latex_options = "repeat_header", font_size = 9) %>%
  add_header_above(c(" " = 2, "Female" = 2, "Male" = 2)) %>%
  row_spec(summ_rows, italic = TRUE) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "top") %>%
  save_kable(file.path(plot.dir, "obs_vs_dapc.tex"))
