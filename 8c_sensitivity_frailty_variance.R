###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Script:   Sensitivity of the baseline (frailty-adjusted) dAPC fit to the
##           assumed frailty variance s2. Total fits (*bs10w) are reused
##           unchanged; only baseline rates depend on s2. Produces the
##           Total-vs-Baseline cohort and centered-period overlays.
##             s2 = 0.50 -> *indiv05.csv,  cohindivnw05.pdf, perplotnew05.pdf
##             s2 = 0.75 -> *indiv075.csv, cohindivnw1.pdf,  perplotnew1.pdf
## Inputs:   <data.dir>/lexis.csv, <data.dir>/multip.csv,
##           cohorteffectbs10w.csv, pereffectbs10w.csv, driftbs10w.csv (Total, from script 02)
## Outputs:  <data.dir> baseline effect/drift/anova CSVs + fits_{suffix}.rds; <plot.dir> the overlays above
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

specs <- list(
  list(s2 = 0.50, suffix = "indiv05",  coh = "cohindivnw05.pdf", per = "perplotnew05.pdf"),
  list(s2 = 0.75, suffix = "indiv075", coh = "cohindivnw1.pdf",  per = "perplotnew1.pdf")
)

## ---- baseline fit at a given s2 (script 4, parameterised) ------------------
multip <- read.table(file.path(data.dir, "multip.csv"), sep = ",", header = TRUE) %>%
  rename(year = Year, age = Age) %>%
  mutate(sex = ifelse(sex == "fem", "Female", "Male"))
lexis  <- read.table(file.path(data.dir, "lexis.csv"), sep = ",", header = TRUE)

fit_baseline <- function(s2, suffix) {
  pc <- lexis %>%
    mutate(year = floor(Year), age = floor(Age)) %>%
    left_join(multip, by = c("country", "sex", "year", "age")) %>%
    mutate(Dxi = Dx * (multip^(-s2)), D = round(Dxi), Y = round(Exp), A = Age, P = Year)
  if (anyNA(pc$D))
    stop(sprintf("NA baseline deaths at s2=%.2f -- check year/age/sex matching.", s2))
  
  countries <- unique(pc$country); sexes <- unique(pc$sex)
  age_list <- per_list <- coh_list <- drift_list <- anova_list <- fits <- list()
  k <- 0L
  for (cc in countries) for (sx in sexes) {
    k <- k + 1L
    dat <- pc %>% filter(country == cc, sex == sx) %>% transmute(A, P, D, Y) %>% as.data.frame()
    m <- apc.fit(data = dat, npar = c(A = 10, P = 10, C = 10),
                 model = "bs", parm = "APC", dr.extr = "Y", ref.p = 1980 + 1/3)
    age_list[[k]]   <- data.frame(m$Age, country = cc, sex = sx, check.names = FALSE)
    per_list[[k]]   <- data.frame(m$Per, country = cc, sex = sx, check.names = FALSE)
    coh_list[[k]]   <- data.frame(m$Coh, country = cc, sex = sx, check.names = FALSE)
    drift_list[[k]] <- data.frame(as.list(m$Drift[1, ]), country = cc, sex = sx, check.names = FALSE)
    anova_list[[k]] <- data.frame(m$Anova, country = cc, sex = sx, check.names = FALSE)
    fits[[paste(cc, sx, sep = "_")]] <- m
  }
  saveRDS(fits, file.path(data.dir, sprintf("fits_%s.rds", suffix)))
  write.table(bind_rows(age_list),   file.path(data.dir, sprintf("ageeffect%s.csv",    suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(per_list),   file.path(data.dir, sprintf("pereffect%s.csv",    suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(coh_list),   file.path(data.dir, sprintf("cohorteffect%s.csv", suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(drift_list), file.path(data.dir, sprintf("drift%s.csv",        suffix)), sep = ",", row.names = FALSE)
  write.table(bind_rows(anova_list), file.path(data.dir, sprintf("myanova%s.csv",      suffix)), sep = ",", row.names = FALSE)
  cat(sprintf("s2=%.2f baseline fitted: %d models -> %s\n", s2, length(fits), suffix))
}

## ---- Total-rates effects (s2-independent; read once) ----------------------
rd <- function(f) read.table(file.path(data.dir, f), sep = ",", header = TRUE) %>%
  distinct() %>% recode_country()
cohTotal   <- rd("cohorteffectbs10w.csv")
perTotal   <- rd("pereffectbs10w.csv")
driftTotal <- rd("driftbs10w.csv") %>% rename(drift = `exp.Est..`)

## ---- plot builders: cohort overlay (A18) and centered period overlay (A19) -
coh_overlay <- function(cohTot, cohBase) {
  d <- bind_rows(cohTot %>% mutate(Rates = "Total"),
                 cohBase %>% mutate(Rates = "Baseline")) %>%
    mutate(Rates = factor(Rates, levels = c("Total", "Baseline")))
  ggplot(d, aes(x = Coh, y = C.RR, color = sexsel, linetype = Rates)) +
    geom_line(size = 0.7) +
    facet_wrap(~ countrysel) +
    scale_x_continuous(
      breaks = seq(1890, 1950, by = 10), limits = c(1886, 1954),
      labels = function(x) paste0(as.integer(x),
                                  "<br><span style='color:grey50;'>", as.integer(x) + 65, "</span>")) +
    scale_y_continuous(limits = c(0.5, 1.4), breaks = seq(0.5, 1.4, by = 0.1)) +  # wider range for higher s2
    geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    labs(x = "Cohort / <span style='color:grey50;'>Year When 65</span>",
         y = "Relative Rate", color = "Sex", linetype = "Rates") +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.x = element_markdown(lineheight = 0.95, margin = margin(t = 3)),
          axis.title.x = element_markdown(size = 12, margin = margin(t = 8)),
          axis.title.y = element_text(size = 12),
          legend.title = element_text(size = 12), legend.position = "bottom",
          panel.spacing.x = unit(1.2, "lines"))
}

make_period_df <- function(per_df, drift_df, rates_label, ref_year = 1980) {
  per_df %>%
    left_join(drift_df, by = c("countrysel", "sexsel")) %>%
    mutate(g_log = log(as.numeric(P.RR)) - (Per - ref_year) * log(drift), Per = floor(Per)) %>%
    group_by(Per, countrysel, sexsel) %>%
    summarise(g_log = mean(g_log, na.rm = TRUE), .groups = "drop") %>%
    mutate(newP.RR = exp(g_log), Rates = rates_label)
}

per_overlay <- function(perTot, driftTot, perBase, driftBase) {
  pereffbs <- bind_rows(
    make_period_df(perTot,  driftTot,  "Total"),
    make_period_df(perBase, driftBase, "Baseline")) %>%
    mutate(Rates = factor(Rates, levels = c("Total", "Baseline")),
           newP.RR = if_else(Per == 1980, 1, newP.RR))
  ## centering: divide each series by its own mean RR (A10/A19 convention)
  sumper <- pereffbs %>% group_by(Rates, countrysel, sexsel) %>%
    summarise(meanRR = mean(newP.RR), .groups = "drop")
  d <- pereffbs %>% left_join(sumper, by = c("Rates","countrysel","sexsel")) %>%
    mutate(RR = newP.RR / meanRR)
  ggplot(d, aes(x = Per, y = RR, color = sexsel, linetype = Rates)) +
    geom_line(size = 0.7) +
    facet_wrap(~ countrysel) +
    labs(x = "Year", y = "Relative Rate", color = "Sex", linetype = "Rates") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10), axis.title = element_text(size = 12),
          legend.title = element_text(size = 12), legend.position = "bottom") +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    scale_x_continuous(breaks = c(seq(1980, 2015, by = 10), 2019), limits = c(1980, 2019)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black")
}

## ---- run each spec: fit baseline, read effects back, plot ------------------
for (s in specs) {
  fit_baseline(s$s2, s$suffix)
  
  cohBase   <- rd(sprintf("cohorteffect%s.csv", s$suffix))
  perBase   <- rd(sprintf("pereffect%s.csv",    s$suffix))
  driftBase <- rd(sprintf("drift%s.csv",        s$suffix)) %>% rename(drift = `exp.Est..`)
  
  ggsave(file.path(plot.dir, s$coh), coh_overlay(cohTotal, cohBase),                         width = 10, height = 12)
  ggsave(file.path(plot.dir, s$per), per_overlay(perTotal, driftTotal, perBase, driftBase),  width = 10, height = 10)
  cat("s2 =", s$s2, "plotted:", s$coh, s$per, "\n")
}
