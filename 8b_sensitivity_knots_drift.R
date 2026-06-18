###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Script:   Sensitivity of the Total-rates dAPC fit to spline knots and drift
##           weighting. Fits two alternative specifications and produces the
##           cohort / detrended-period / age panels for each.
##             (a) 12 knots, exposure-weighted drift -> *bs12w.csv, cohlpt12w.pdf, perplotnew12w.pdf, ageplot12w.pdf
##             (b) 10 knots, non-weighted drift       -> *bs10n.csv, cohlpt10n.pdf, perplotnew10n.pdf, ageplot10n.pdf
##           Point estimates + model-based CI columns retained (as in script 02); no Monte Carlo CI.
## Input:    <data.dir>/lexis.csv
## Outputs:  <data.dir> effect/drift/anova CSVs + fits_bs{suffix}.rds; <plot.dir> the panels above
## Depends:  Epi, dplyr, tidyr, ggplot2, ggtext
###############################################################################
rm(list = ls())

library(Epi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)

data.dir   <- #"your/data/directory"   
plot.dir   <- #"your/plots/directory" 
  

cntry_map <- c(CHE = "Switzerland", DNK = "Denmark",FIN = "Finland",FRATNP = "France",GBRTENW = "England and Wales",  
               ITA = "Italy",NLD = "Netherlands", NOR = "Norway",SWE = "Sweden")

recode_country <- function(df) {
  if (!"countrysel" %in% names(df) && "country" %in% names(df))
    df <- dplyr::rename(df, countrysel = country)
  if (!"sexsel"     %in% names(df) && "sex"     %in% names(df))
    df <- dplyr::rename(df, sexsel = sex)
  df %>% mutate(countrysel = recode(countrysel, !!!as.list(cntry_map),
                                    .default = countrysel))
}

lexis     <- read.table(file.path(data.dir, "lexis.csv"), sep = ",", header = TRUE)
countries <- unique(lexis$country)
sexes     <- unique(lexis$sex)

## ---- the two sensitivity specifications ------------------------------------
## knots = npar per dimension; drextr = drift-weighting passed to apc.fit
specs <- list(
  list(suffix = "12w", knots = 12, drextr = "Y"),    # 12 knots, exposure-weighted
  list(suffix = "10n", knots = 10, drextr = "")       # 10 knots, non-weighted
)

## ---- fit one spec: one model per country x sex, accumulate, write CSVs ------
fit_spec <- function(knots, drextr, suffix) {
  age_list <- per_list <- coh_list <- drift_list <- anova_list <- fits <- list()
  k <- 0L
  for (cc in countries) {
    for (sx in sexes) {
      k <- k + 1L
      dat <- lexis %>%
        filter(country == cc, sex == sx) %>%
        transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>%
        as.data.frame()
      m <- apc.fit(data = dat,
                   npar = c(A = knots, P = knots, C = knots),
                   model = "bs", parm = "APC",
                   dr.extr = drextr, ref.p = 1980 + 1/3)
      fits[[paste(cc, sx, sep = "_")]] <- m
      age_list[[k]]   <- data.frame(m$Age, country = cc, sex = sx, check.names = FALSE)
      per_list[[k]]   <- data.frame(m$Per, country = cc, sex = sx, check.names = FALSE)
      coh_list[[k]]   <- data.frame(m$Coh, country = cc, sex = sx, check.names = FALSE)
      drift_list[[k]] <- data.frame(as.list(m$Drift[1, ]), country = cc, sex = sx, check.names = FALSE)
      anova_list[[k]] <- data.frame(m$Anova, country = cc, sex = sx, check.names = FALSE)
    }
  }
  saveRDS(fits, file.path(data.dir, sprintf("fits_bs%s.rds", suffix)))
  write.table(bind_rows(age_list),   file.path(data.dir, sprintf("ageeffectbs%s.csv",    suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(per_list),   file.path(data.dir, sprintf("pereffectbs%s.csv",    suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(coh_list),   file.path(data.dir, sprintf("cohorteffectbs%s.csv", suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(drift_list), file.path(data.dir, sprintf("driftbs%s.csv",        suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(anova_list), file.path(data.dir, sprintf("myanovabs%s.csv",      suffix)), sep = ",", row.names = FALSE)
  cat("spec", suffix, "fitted:", length(fits), "models\n")
}

## ---- reusable plot builders (Total rates only; colour = sex) ---------------
coh_plot <- function(coh) {
  ggplot(coh, aes(x = Coh, y = C.RR, color = sexsel)) +
    geom_line(size = 0.7) +
    facet_wrap(~ countrysel) +
    scale_x_continuous(
      breaks = seq(1890, 1950, by = 10), limits = c(1886, 1954),
      labels = function(x) paste0(as.integer(x),
                                  "<br><span style='color:grey50;'>", as.integer(x) + 65, "</span>")) +
    scale_y_continuous(limits = c(0.5, 1.15), breaks = seq(0.5, 1.1, by = 0.1)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    labs(x = "Cohort / <span style='color:grey50;'>Year When 65</span>",
         y = "Relative Rate", color = "Sex") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.x = element_markdown(lineheight = 0.95, margin = margin(t = 3)),
          axis.title.x = element_markdown(size = 12, margin = margin(t = 8)),
          axis.title.y = element_text(size = 12),
          legend.title = element_text(size = 12), legend.position = "bottom",
          panel.spacing.x = unit(1.2, "lines"))
}

make_period_df <- function(per_df, drift_df, ref_year = 1980) {
  per_df %>%
    left_join(drift_df, by = c("countrysel", "sexsel")) %>%
    mutate(g_log = log(as.numeric(P.RR)) - (Per - ref_year) * log(drift),
           Per = floor(Per)) %>%
    group_by(Per, countrysel, sexsel) %>%
    summarise(g_log = mean(g_log, na.rm = TRUE), .groups = "drop") %>%
    mutate(newP.RR = exp(g_log),
           newP.RR = if_else(Per == ref_year, 1, newP.RR))
}

per_plot <- function(per_df, drift_df) {
  d <- make_period_df(per_df, drift_df)
  ggplot(d, aes(x = Per, y = newP.RR, color = sexsel)) +
    geom_line(size = 0.7) +
    facet_wrap(~ countrysel) +
    labs(x = "Year", y = "Relative Rate", color = "Sex") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10), axis.title = element_text(size = 12),
          legend.title = element_text(size = 12), legend.position = "bottom") +
    scale_x_continuous(breaks = c(seq(1980, 2015, by = 10), 2019), limits = c(1980, 2019)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black")
}

age_plot <- function(age) {
  ggplot(age, aes(x = Age, y = log(Rate), color = sexsel)) +
    geom_line(size = 0.7) +
    facet_wrap(~ countrysel) +
    labs(x = "Age", y = "Mortality rate (log-scale)", color = "Sex") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10), axis.title = element_text(size = 12),
          legend.title = element_text(size = 12)) +
    scale_x_continuous(breaks = seq(65, 95, by = 10), limits = c(65, 95))
}

## ---- run each spec: fit, then read effects back and plot -------------------
rd <- function(f) read.table(file.path(data.dir, f), sep = ",", header = TRUE) %>%
  distinct() %>% recode_country()

for (s in specs) {
  fit_spec(s$knots, s$drextr, s$suffix)
  
  coh   <- rd(sprintf("cohorteffectbs%s.csv", s$suffix))
  per   <- rd(sprintf("pereffectbs%s.csv",    s$suffix))
  age   <- rd(sprintf("ageeffectbs%s.csv",    s$suffix))
  drift <- rd(sprintf("driftbs%s.csv",        s$suffix)) %>% rename(drift = `exp.Est..`)
  
  ggsave(file.path(plot.dir, sprintf("cohlpt%s.pdf",     s$suffix)), coh_plot(coh),        width = 10, height = 12)
  ggsave(file.path(plot.dir, sprintf("perplotnew%s.pdf", s$suffix)), per_plot(per, drift), width = 10, height = 10)
  ggsave(file.path(plot.dir, sprintf("ageplot%s.pdf",    s$suffix)), age_plot(age),        width = 10, height = 10)
  cat("spec", s$suffix, "plotted\n")
}
