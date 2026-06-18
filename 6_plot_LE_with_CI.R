###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Draw Figures 1-3 from the Monte Carlo CI table, with 95%
##           ribbons on every model-derived series. Central lines and bands both
##           come from LE_series_ci.csv (05), so they are mutually consistent.
## Inputs:   <models.dir>/LE_series_ci.csv   (from 05; country/sex are HMD codes)
##           <out.dir>/lexis.csv           (observed mx, for the dots)
## Outputs:  <plot.dir>/LE_impr_absolute.pdf, LE_impr_decomp.pdf, LE_impr_and_gap_absolute.pdf
## Depends:  dplyr, tidyr, ggplot2, patchwork, slider
###############################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(slider)

data.dir   <- #"your/data/directory"  
models.dir <- #"your/output/directory" 
plot.dir   <- #"your/plots/directory" 

ages  <- 65:94
years <- 1980:2019
yr_lo <- 1992; yr_hi <- 2018             # plotted window (year > yr_lo & < yr_hi)

cntry_map <- c(CHE = "Switzerland", DNK = "Denmark", FIN = "Finland", FRATNP = "France",
               GBRTENW = "England & Wales", ITA = "Italy", NLD = "Netherlands",
               NOR = "Norway", SWE = "Sweden")

## same partial-LE engine as 05, so observed dots use the same life-table closure
pLE_matrix <- function(M) {
  K <- nrow(M)
  qx <- M / (1 + 0.5 * M); qx[K, ] <- 1
  px <- 1 - qx
  lx <- matrix(1, K, ncol(M))
  if (K > 1) for (k in 2:K) lx[k, ] <- lx[k - 1, ] * px[k - 1, ]
  dx <- matrix(0, K, ncol(M))
  if (K > 1) dx[1:(K - 1), ] <- lx[1:(K - 1), ] - lx[2:K, ]
  dx[K, ] <- lx[K, ]
  Lx <- matrix(0, K, ncol(M))
  if (K > 1) Lx[1:(K - 1), ] <- lx[2:K, ] + 0.5 * dx[1:(K - 1), ]
  Lx[K, ] <- lx[K, ] / M[K, ]
  colSums(Lx)
}

## ---- read the CI table -----------------------------------------------------
ci <- read.csv(file.path(models.dir, "LE_series_ci.csv")) %>%
  mutate(country = recode(country, !!!as.list(cntry_map)))

## ---- observed Total improvement (for the Figure 1 dots) --------------------
obs_impr <- read.table(file.path(models.dir, "lexis.csv"), sep = ",", header = TRUE) %>%
  mutate(year = floor(Year), age = floor(Age)) %>%
  filter(age %in% ages, year %in% years) %>%
  group_by(country, sex, year, age) %>%
  summarise(Dx = sum(Dx), Exp = sum(Exp), .groups = "drop") %>%   # pool the two triangles
  mutate(mx = Dx / Exp) %>%
  arrange(country, sex, year, age) %>%
  group_by(country, sex, year) %>%
  summarise(LE = pLE_matrix(matrix(mx, ncol = 1))[1], .groups = "drop") %>%
  arrange(country, sex, year) %>%
  group_by(country, sex) %>%
  mutate(d1 = LE - lag(LE),
         impr = slide_dbl(d1, mean, .before = 9, .complete = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(impr), year > yr_lo, year < yr_hi) %>%
  mutate(country = recode(country, !!!as.list(cntry_map)))


## =============================================================================
## FIGURE 1: component improvements (Total) with ribbons + observed dots
## =============================================================================
f1 <- ci %>%
  filter(rates == "Total", series %in% c("dA", "dAC", "dAP", "dAPC"),
         year > yr_lo, year < yr_hi) %>%
  mutate(series = factor(series, levels = c("dAPC", "dA", "dAC", "dAP")))

pal1 <- c(dAPC = "#1b9e77", dA = "#d95f02", dAC = "#7570b3", dAP = "#e7298a")
f1_rng <- range(c(f1$lower, f1$upper, obs_impr$impr), na.rm = TRUE)

fig1_panel <- function(sx) {
  ggplot(f1 %>% filter(sex == sx),
         aes(year, estimate, color = series, fill = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, colour = NA) +
    geom_line(aes(y = lower), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(aes(y = upper), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(linewidth = 0.7) +
    geom_point(data = obs_impr %>% filter(sex == sx),
               aes(year, impr), inherit.aes = FALSE, size = 0.5, colour = "grey20") +
    facet_wrap(~ country) +
    scale_color_manual(values = pal1) + scale_fill_manual(values = pal1) +
    scale_x_continuous(breaks = c(seq(1990, 2010, by = 10), 2019), limits = c(1990, 2019)) +
    coord_cartesian(ylim = f1_rng) +
    labs(x = "Year", y = "Change (in Years)", color = NULL, fill = NULL,
         subtitle = sx) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.caption  = element_text(size = 9, hjust = 0),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10), axis.title = element_text(size = 12),
          legend.position = "bottom")
}
fig1 <- (fig1_panel("Female") + fig1_panel("Male")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(plot.dir, "LE_impr_absolute.pdf"), fig1, width = 17, height = 10)




## =============================================================================
## FIGURE 2: total improvement (dAPC) decomposed into its additive
## Shapley parts — drift (dA), non-linear period, non-linear cohort — with the
## observed improvement as dots. dA + period + cohort = dAPC at every year.
## =============================================================================
f2b <- ci %>%
  filter(rates == "Total",
         series %in% c("dAPC", "dA", "shap_period", "shap_cohort"),
         year > yr_lo, year < yr_hi) %>%
  mutate(series = recode(series,
                         dAPC        = "Total",
                         dA          = "Drift",
                         shap_period = "Non-linear period",
                         shap_cohort = "Non-linear cohort"),
         series = factor(series,
                         levels = c("Total", "Drift",
                                    "Non-linear period", "Non-linear cohort")))

## consistency: drift + period + cohort must equal the total at every year
chk <- ci %>%
  filter(rates == "Total",
         series %in% c("dA", "shap_period", "shap_cohort", "dAPC"),
         year > yr_lo, year < yr_hi) %>%
  dplyr::select(country, sex, year, series, estimate) %>%
  tidyr::pivot_wider(names_from = series, values_from = estimate) %>%
  mutate(resid = dA + shap_period + shap_cohort - dAPC)
stopifnot(max(abs(chk$resid), na.rm = TRUE) < 1e-6)

pal2b <- c("Total"       = "#111111",
           "Drift"         = "#d95f02",
           "Non-linear period"  = "#e7298a",
           "Non-linear cohort"  = "#7570b3")

f2b_rng <- range(c(f2b$lower, f2b$upper, obs_impr$impr), na.rm = TRUE)

fig2b_panel <- function(sx) {
  ggplot(f2b %>% filter(sex == sx),
         aes(year, estimate, color = series, fill = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, colour = NA) +
    geom_line(aes(y = lower), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(aes(y = upper), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(linewidth = 0.7) +
    geom_point(data = obs_impr %>% filter(sex == sx),
               aes(year, impr), inherit.aes = FALSE, size = 0.5, colour = "grey20") +
    facet_wrap(~ country) +
    scale_color_manual(values = pal2b) + scale_fill_manual(values = pal2b) +
    scale_x_continuous(breaks = c(seq(1990, 2010, by = 10), 2019), limits = c(1990, 2019)) +
    coord_cartesian(ylim = f2b_rng) +
    labs(x = "Year", y = "Change / Contribution (in Years)",
         color = NULL, fill = NULL, subtitle = sx) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.caption  = element_text(size = 9, hjust = 0),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10), axis.title = element_text(size = 12),
          legend.position = "bottom")
}

fig2b <- (fig2b_panel("Female") + fig2b_panel("Male")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(plot.dir, "LE_impr_decomp.pdf"), fig2b, width = 17, height = 10)


## =============================================================================
## FIGURE 3: improvements (Total vs Baseline) + gap, with ribbons
## =============================================================================
le_panel <- ci %>%
  filter(series == "dAPC", rates %in% c("Total", "Baseline")) %>%
  transmute(country, sex, year, Line = rates, estimate, lower, upper,
            Panel = "LE improvements")

gap_panel <- ci %>%
  filter(series == "gap_dAPC") %>%
  transmute(country, sex, year,
            Line = "Gap",
            estimate, lower, upper, Panel = "Gap")

plot_df <- bind_rows(le_panel, gap_panel) %>%
  filter(year > yr_lo, year < yr_hi) %>%
  mutate(Line  = factor(Line,  levels = c("Total", "Baseline", "Gap")),
         Panel = factor(Panel, levels = c("LE improvements", "Gap")))

pal2  <- c(Total = "#1b9e77", Baseline = "#d95f02", Gap = "#7570b3")
lims2 <- list(
  `LE improvements` = range(c(plot_df$lower[plot_df$Panel == "LE improvements"],
                              plot_df$upper[plot_df$Panel == "LE improvements"]), na.rm = TRUE),
  `Gap`             = range(c(plot_df$lower[plot_df$Panel == "Gap"],
                              plot_df$upper[plot_df$Panel == "Gap"]), na.rm = TRUE)
)

make_sex_panel_plot <- function(sex_label, panel_label) {
  ggplot(plot_df %>% filter(sex == sex_label, Panel == panel_label),
         aes(year, estimate, color = Line, fill = Line)) +
    {if (panel_label == "Gap")
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50")} +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
    geom_line(aes(y = lower), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(aes(y = upper), linewidth = 0.25, linetype = "22", alpha = 0.7) +
    geom_line(linewidth = 0.8) +
    facet_grid(cols = vars(country)) +
    scale_color_manual(values = pal2, name = "") +
    scale_fill_manual(values = pal2, name = "") +
    coord_cartesian(ylim = lims2[[panel_label]]) +
    scale_x_continuous(breaks = c(seq(1990, 2010, by = 10), 2019), limits = c(1990, 2019)) +
    labs(x = "Year", y = "Change in Years", subtitle = panel_label) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 9), axis.title = element_text(size = 11),
          panel.spacing = unit(0.6, "lines"), legend.position = "bottom")
}

make_sex_plot <- function(sex_label) {
  p_le  <- make_sex_panel_plot(sex_label, "LE improvements") +
    ggtitle(sex_label) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  p_gap <- make_sex_panel_plot(sex_label, "Gap")
  (p_le / p_gap) + plot_layout(heights = c(1, 1))
}

fig2 <- wrap_plots(make_sex_plot("Female"), make_sex_plot("Male"),
                   ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(plot.dir, "LE_impr_and_gap_absolute.pdf"), fig2, width = 17, height = 11)