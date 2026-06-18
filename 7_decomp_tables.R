###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Supplementary tables of the Shapley period decomposition, ONE TABLE
##           PER SEX, all values rounded to two decimals (estimate and CI).
##           Reads the numeric long file (not the pre-formatted wide one) so the
##           2-dp rounding is done once, cleanly.
## Input:    <models.dir>/LE_change_decomp.csv   (from 05; long: estimate/lower/upper)
## Outputs:  <plot.dir>/decomp_table_female.tex, <plot.dir>/decomp_table_male.tex
## Depends:  dplyr, tidyr, kableExtra
###############################################################################
library(dplyr); 
library(tidyr); 
library(kableExtra)

  models.dir <- #"your/data/directory" 
  plot.dir   <- #"your/plots/directory" 
  

cntry_map <- c(CHE = "Switzerland", DNK = "Denmark", FIN = "Finland", FRATNP = "France",
               GBRTENW = "England & Wales", ITA = "Italy", NLD = "Netherlands",
               NOR = "Norway", SWE = "Sweden")
period_lv <- c("1983-1993", "1993-2003", "2003-2013", "2013-2017")

## format "est (lo, hi)" at 2 dp; %.2f keeps trailing zeros (0.10, not 0.1)
fmt2 <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)

dec <- read.csv(file.path(models.dir, "LE_change_decomp.csv")) %>%
  filter(component %in% c("drift", "period", "cohort", "total")) %>%
  mutate(country   = recode(country, !!!as.list(cntry_map)),
         period    = factor(period, levels = period_lv),
         component = recode(component,
                            drift  = "Drift",
                            period = "NL period",
                            cohort = "NL cohort",
                            total  = "Total"),
         component = factor(component,
                            levels = c("Drift", "NL period", "NL cohort", "Total")),
         cell      = fmt2(estimate, lower, upper)) %>%
  arrange(country, sex, period, component)

## one wide table per sex, then write each to its own .tex
make_sex_table <- function(sx, fname) {
  w <- dec %>%
    filter(sex == sx) %>%
    dplyr::select(country, period, component, cell) %>%
    pivot_wider(names_from = component, values_from = cell) %>%
    arrange(country, period) %>%
    mutate(period = gsub("-", "\u2013", period)) %>%        # en-dash
    dplyr::select(Country = country, Period = period,
                  Drift, `NL period`, `NL cohort`, Total)
  
  kbl(w, format = "latex", booktabs = TRUE, escape = TRUE, linesep = "",
      align = c("l", "l", "r", "r", "r", "r"),
      caption = sprintf(paste("Shapley decomposition of the average annual change",
                              "in partial life expectancy (ages 65--94), %s, in",
                              "years per year with 95\\%% Monte Carlo intervals."),
                        tolower(sx)),
      label = sprintf("decomp_%s", tolower(sx))) %>%
    kable_styling(latex_options = c("repeat_header"), font_size = 9) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "top") %>%
    save_kable(file.path(plot.dir, fname))
}

make_sex_table("Female", "decomp_table_female.tex")
make_sex_table("Male",   "decomp_table_male.tex")