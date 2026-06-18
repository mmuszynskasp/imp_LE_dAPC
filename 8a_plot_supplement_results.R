###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Script:   Supplementary exhibits computed on the 1980-2019 window, plus the
##           1960-vs-1980 cohort-effect overlay and the drift-only restricted
##           e-dagger / entropy supplement (Hypothesis 1).
## Inputs:   <data.dir> : effect/drift/anova CSVs (*bs10w, *indiv10w; from 02/04),
##                        cohorteffectbs1960.csv / pereffectbs1960.csv / driftbs1960.csv (from 02_1960),
##                        dA_edag_H.csv (from 05), pc_data.csv & pc_data_si.csv (observed / standard rates)
## Outputs:  <plot.dir> : anovaplot.pdf, cohlpt_both.pdf, perplotnew.pdf, perplotnew2.pdf,
##                        ageplot.pdf, exAPC.pdf, exAPCi.pdf, eximpr10w05s.pdf,
##                        cohort_effect_1960_vs_1980.pdf, edag_H_dA.pdf,
##                        + LaTeX anova/drift tables printed to console (xtable)
## Depends:  tidyr, dplyr, ggplot2, ggpubr, patchwork, mgcv, MortalityLaws, xtable, ggtext
###############################################################################
rm(list = ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(mgcv)
library(MortalityLaws)
library(xtable)
library(ggtext)

data.dir   <- #"your/data/directory"  
plot.dir   <- #"your/plots/directory" 
  

ourcohorts <- 1885:1954
years <- 1980:2019
ages  <- 65:94

#########read in data from both models (1980-2019 main window) #################

cohTotal <- read.table(file=file.path(data.dir, "cohorteffectbs10w.csv"),sep=",",header=TRUE) %>%
  distinct()
cohIndiv <- read.table(file=file.path(data.dir, "cohorteffectindiv10w.csv"),sep=",",header=TRUE)%>%
  distinct()
perTotal <- read.table(file=file.path(data.dir, "pereffectbs10w.csv"),sep=",",header=TRUE)%>%
  distinct()
perIndiv <- read.table(file=file.path(data.dir, "pereffectindiv10w.csv"),sep=",",header=TRUE)%>%
  distinct()
ageTotal <- read.table(file=file.path(data.dir, "ageeffectbs10w.csv"),sep=",",header=TRUE)%>%
  distinct()
ageIndiv <- read.table(file=file.path(data.dir, "ageeffectindiv10w.csv"),sep=",",header=TRUE)%>%
  distinct()
driftTotal <- read.table(file=file.path(data.dir, "driftbs10w.csv"),sep=",", header=TRUE) %>%
  distinct()
driftIndiv <- read.table(file=file.path(data.dir, "driftindiv10w.csv"),sep=",", header=TRUE)%>%
  distinct()

driftTotal <- driftTotal %>% rename(drift = `exp.Est..`)
driftIndiv <- driftIndiv %>% rename(drift = `exp.Est..`)

anovaTotal <- read.table(file=file.path(data.dir, "myanovabs10w.csv"),sep=",",header=TRUE) %>% 
  rename("countrysel"="country") #for replacement of country names in the next step
anovaIndiv <- read.table(file=file.path(data.dir, "myanovaindiv10w.csv"),sep=",",header=TRUE)%>% 
  rename("countrysel"="country")#for replacement of country names in the next step

###recode country names for the plots
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
objs <- c("cohTotal","cohIndiv","perTotal","perIndiv","ageTotal","ageIndiv","driftTotal","driftIndiv","anovaTotal","anovaIndiv")

updated <- lapply(mget(objs, inherits = TRUE), recode_country)
list2env(updated, .GlobalEnv)

###############################################################################
## 1. Contribution of effects (deviance-based %)  -> anovaplot + tables
###############################################################################
totest <- anovaTotal %>%
  select(countrysel, sexsel, Model, `Mod..dev.`) %>%
  mutate(Type = "Total") %>%
  add_row(anovaIndiv %>%
            select(countrysel, sexsel, Model, `Mod..dev.`) %>%
            mutate(Type = "Baseline")) %>%
  pivot_wider(names_from = Model, values_from = `Mod..dev.`, values_fn = dplyr::first) %>%
  mutate(country = countrysel,
         reddrift = Age - `Age-drift`,
         redP     = `Age-drift` - `Age-Period`,
         redC     = `Age-Period` - `Age-Period-Cohort`,
         total    = Age - `Age-Period-Cohort`,
         Drift     = 100 * reddrift / total,
         Period    = 100 * redP / total,
         Cohort    = 100 * redC / total,
         Total_exp = 100 * total / Age) %>%
  select(country, sex = sexsel, Cohort, Period, Total_exp, Type) %>%   # rename back to `sex` for downstream
  pivot_longer(Cohort:Total_exp, names_to = "Effect", values_to = "Contribution")

myaverages <- totest %>%
  group_by(sex,Type,Effect) %>%
  summarise(Contribution=mean(Contribution)) %>%
  mutate(country="Average")

type_order <- c("Baseline", "Total")
totest1 <- totest %>%
  add_row(myaverages) %>%
  mutate(Type = factor(Type, levels = type_order))

#separate the effects
pc_data <- totest1 %>%
  filter(Effect %in% c("Period", "Cohort"))

drift_data <- totest1 %>%
  filter(Effect == "Drift") %>%
  group_by(country, sex, Type) %>%
  summarise(Drift = first(Contribution), .groups = "drop")


pc_top <- pc_data %>%
  group_by(country, sex, Type) %>%
  summarise(top = sum(Contribution, na.rm = TRUE), .groups = "drop")

labels_df <- pc_top %>%
  left_join(drift_data, by = c("country", "sex", "Type")) %>%
  mutate(Drift = ifelse(is.na(Drift), 100 - top, Drift))

country_levels <- levels(factor(totest1$country))
top_country <- country_levels[1]
type_levels <- levels(totest1$Type)             
top_type <- type_levels[length(type_levels)] 
drift_header <- totest1 %>%
  distinct(sex) %>%
  mutate(country = top_country, Type = top_type, x = 32, label = "Drift")


anovaplot <- ggplot() +
  geom_col(data = pc_data,aes(y = Type, x = Contribution, fill = Effect), position = "stack") +
  geom_text(data = labels_df, aes(y = Type, x = 32, label = sprintf("%.0f%%", Drift)),size = 3.2, hjust = 0.5, na.rm = TRUE) +
  geom_text(data = drift_header,aes(x = x, y = Type, label = label),inherit.aes = FALSE,vjust = -1.1, size = 3.8, fontface = "bold") +
  facet_grid(country ~ sex, scales = "free_y", space = "free_y") +
  labs(y = NULL, x = "Contribution (%)", fill = "Effect") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.6))) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(plot.margin = margin(14, 12, 6, 6),
    strip.text   = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 11),
    axis.title   = element_text(size = 12),
    legend.title = element_text(size = 12))


ggsave(anovaplot,filename=file.path(plot.dir, "anovaplot.pdf"),width=11,height=14)


### table (Female rows 1-9, Male rows 10-18)
tab_anova <- totest %>% filter(sex=="Female") %>%
  pivot_wider(names_from = c(Effect,Type), values_from = Contribution, names_glue = "{Effect}|{Type}") %>%
  arrange(country) %>%
  mutate(`Drift|Total`=100-`Cohort|Total`-`Period|Total`,`Drift|Baseline`=100-`Cohort|Baseline`-`Period|Baseline`) %>%
  select(country,`Cohort|Total`,`Cohort|Baseline`,`Period|Total`,`Period|Baseline`,`Drift|Total`,`Drift|Baseline`,
         `Total_exp|Total`,`Total_exp|Baseline`) %>%
  add_row(totest %>% filter(sex=="Male") %>%
    pivot_wider(names_from = c(Effect,Type), values_from = Contribution, names_glue = "{Effect}|{Type}") %>%
    mutate(`Drift|Total`=100-`Cohort|Total`-`Period|Total`,`Drift|Baseline`=100-`Cohort|Baseline`-`Period|Baseline`) %>%
    select(country,`Cohort|Total`,`Cohort|Baseline`,`Period|Total`,`Period|Baseline`,`Drift|Total`,`Drift|Baseline`,
             `Total_exp|Total`,`Total_exp|Baseline`))

print(xtable(tab_anova, digits=1), include.rownames = FALSE)

tab_anova_averF <- as.data.frame(t(colMeans(tab_anova[1:9,-1])))
print(xtable(tab_anova_averF, digits = 1), include.rownames = FALSE)

tab_anova_averM <- as.data.frame(t(colMeans(tab_anova[10:18,-1])))
print(xtable(tab_anova_averM, digits=1), include.rownames = FALSE)

###############################################################################
## 2. Cohort effect (Total vs Baseline, 1980-2019)  -> cohlpt_both.pdf
###############################################################################
cohboth <- cohTotal %>%
  mutate(Rates="Total") %>%
  add_row(cohIndiv %>%
  mutate(Rates="Baseline")) %>%
  mutate(Rates=factor(Rates, levels=c("Total","Baseline")))

x_breaks <- c(seq(1890, 1950, by = 10))

cohplotI <- ggplot(cohboth, aes(x = Coh, y = C.RR, color = sexsel, linetype = Rates)) +
  geom_line(size = 0.7) +
  facet_wrap(~ countrysel) +
  scale_x_continuous(
    breaks = x_breaks,
    limits = c(1886, 1954),
    labels = function(x) {
      paste0(as.integer(x),"<br><span style='color:grey50;'>", as.integer(x) + 65, 
             "</span>")}) +
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
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    panel.spacing.x = unit(1.2, "lines"))

ggsave(cohplotI,filename=file.path(plot.dir, "cohlpt_both.pdf"),width=10,height=12)


###############################################################################
## 3. Period effect (drift removed, 1980-2019)  -> perplotnew.pdf / perplotnew2.pdf
###############################################################################
library(dplyr)

make_period_df <- function(per_df, drift_df, rates_label, ref_year = 1980) {
  per_df %>%
     left_join(drift_df, by = c("countrysel", "sexsel")) %>%
    mutate(g_log = log(as.numeric(P.RR)) - (Per - ref_year) * log(drift)) %>%
    mutate(Per = floor(Per)) %>%
    group_by(Per, countrysel, sexsel) %>%
    summarise(g_log = mean(g_log, na.rm = TRUE),
              P.RR  = mean(as.numeric(P.RR), na.rm = TRUE), .groups = "drop") %>%
    mutate(newP.RR = exp(g_log), Rates  = rates_label)}

pereffbs <- bind_rows(
  make_period_df(perTotal, driftTotal, "Total", ref_year = 1980),
  make_period_df(perIndiv, driftIndiv, "Baseline", ref_year = 1980)) %>%
  mutate(Rates = factor(Rates, levels = c("Total", "Baseline")),
    # optional: anchor the plotted RR at 1980 = 1 (presentation only)
    newP.RR = if_else(Per == 1980, 1, newP.RR))


####plot without the drift, annual RR
perplotnew <- ggplot(pereffbs%>%
                       filter(Rates=="Total"), aes(x = Per, y = newP.RR, color = sexsel, linetype=Rates)) +
  geom_line(size=0.7)+
  facet_wrap(~ countrysel) +
  labs(
    x = "Year",
    y = "Relative Rate",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12))+
  guides(linetype = guide_legend(override.aes = list(color = "black")))+
  scale_x_continuous(breaks = c(seq(1980, 2015, by = 10), 2019), limits = c(1980, 2019))+
  theme(legend.position = "bottom")+ 
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

ggsave(perplotnew,filename=file.path(plot.dir, "perplotnew.pdf"),width=10,height=10)



sumper <- pereffbs %>%
  group_by(Rates,countrysel,sexsel) %>%
  summarise(meanRR=mean(newP.RR))


sumper <- pereffbs %>%
  left_join(sumper) %>%
  mutate(RR=newP.RR/meanRR)


perplotnew2 <- ggplot(sumper , aes(x = Per, y = RR, color = sexsel, linetype=Rates)) +
  geom_line(size=0.7)+
  facet_wrap(~ countrysel) +
  labs(
    x = "Year",
    y = "Relative Rate",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12))+
  guides(linetype = guide_legend(override.aes = list(color = "black")))+
  scale_x_continuous(breaks = c(seq(1980, 2015, by = 10), 2019), limits = c(1980, 2019))+
  theme(legend.position = "bottom")+ 
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

ggsave(perplotnew2,filename=file.path(plot.dir, "perplotnew2.pdf"),width=10,height=10)


###############################################################################
## 4. Drift table  (drift is the ANNUAL rate ratio - printed raw, no ^0.5)
###############################################################################
drift <- driftTotal %>%
  mutate(Rates = "Total") %>%
  add_row(driftIndiv %>% mutate(Rates = "Baseline")) %>%
  select(countrysel, sexsel, drift, Rates)          # <- drop X2.5. / X97.5.

tabledrift <- drift %>%
  filter(sexsel == "Female") %>%
  select(-sexsel) %>%
  pivot_wider(names_from = Rates, values_from = drift) %>%
  left_join(drift %>%
              filter(sexsel == "Male") %>%
              select(-sexsel) %>%
              mutate(Rates = ifelse(Rates == "Total", "Totalm", "Baselinem")) %>%
              pivot_wider(names_from = Rates, values_from = drift),
            by = "countrysel") %>%
  arrange(countrysel)

library(xtable)
print(xtable(tabledrift, digits=3), include.rownames = FALSE)

averdrift <- colMeans(tabledrift[,-1])
print(xtable(t(averdrift), digits=3), include.rownames = FALSE)

###############################################################################
## 5. Age effect  -> ageplot.pdf
###############################################################################
ageeffbs <- ageTotal  %>%
  distinct() %>%
  mutate(Type="Total") %>%
  add_row(ageIndiv  %>%
            distinct() %>%
            mutate(Type="Baseline"))%>%
  mutate(Type = factor(Type, levels = c("Total", "Baseline")))



ageplot <- ggplot(ageeffbs, aes(x = Age, y = log(Rate), color = sexsel, linetype=Type)) +
  geom_line(size=0.7) +
  facet_wrap(~ countrysel) +
  labs(
    x = "Age",
    y = "Mortality rate (log-scale)",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12))+
  scale_x_continuous(breaks = c(seq(65, 95, by = 10)), limits = c(65, 95))


ggsave(ageplot,filename=file.path(plot.dir, "ageplot.pdf"),width=10,height=10)

###############################################################################
## 6. LE reconstruction from the effects  -> LE plots
###############################################################################
#new rates taking into account drift, non-linear period and age only 
newexTotal <- perTotal %>%
  mutate(Per=floor(Per)) %>%
  group_by(Per,countrysel,sexsel) %>%       ####necessary because Lexis triangles
  summarise(P.RR=mean(as.numeric(P.RR)))%>%
  inner_join(ageTotal%>%
               mutate(Age=floor(as.numeric(Age))) %>%
               group_by(Age,countrysel,sexsel) %>%       ####necessary because Lexis triangles
               summarise(Rate=mean(as.numeric(Rate))), by = c("countrysel", "sexsel")) %>%
  mutate(Coh=Per-Age) %>%
  left_join(cohTotal%>%
               mutate(Coh=floor(as.numeric(Coh))) %>%
               group_by(Coh,countrysel,sexsel) %>%       ####necessary because Lexis triangles
               summarise(C.RR=mean(as.numeric(C.RR)))) %>%
  left_join(driftTotal) %>%
  mutate(dAP=Rate*P.RR,
         dAC=Rate*C.RR,
         dA=Rate*(drift^(Per-1980)),
         dAPC=Rate*C.RR*P.RR) %>%
  group_by(countrysel,sexsel,Per) %>%
  mutate(exdAP=LifeTable(x=Age,mx=dAP,ax=0.5)$lt$ex,
         exdAC=LifeTable(x=Age,mx=dAC,ax=0.5)$lt$ex,
         exdA=LifeTable(x=Age,mx=dA,ax=0.5)$lt$ex,
         exdAPC=LifeTable(x=Age,mx=dAPC,ax=0.5)$lt$ex) %>%
  ungroup() %>%
  filter(Age==65) 


newexIndiv <- perIndiv %>%
  mutate(Per=floor(Per)) %>%
  group_by(Per,countrysel,sexsel) %>%       ####necessary because Lexis triangles
  summarise(P.RR=mean(as.numeric(P.RR)))%>%
  inner_join(ageIndiv%>%
               mutate(Age=floor(as.numeric(Age))) %>%
               group_by(Age,countrysel,sexsel) %>%       ####necessary because Lexis triangles
               summarise(Rate=mean(as.numeric(Rate))), by = c("countrysel", "sexsel")) %>%
  mutate(Coh=Per-Age) %>%
  left_join(cohIndiv%>%
              mutate(Coh=floor(as.numeric(Coh))) %>%
              group_by(Coh,countrysel,sexsel) %>%       ####necessary because Lexis triangles
              summarise(C.RR=mean(as.numeric(C.RR)))) %>%
  left_join(driftIndiv) %>%
  mutate(dAP=Rate*P.RR,
         dAC=Rate*C.RR,
         dA=Rate*(drift^(Per-1980)),
         dAPC=Rate*C.RR*P.RR) %>%
  group_by(countrysel,sexsel,Per) %>%
  mutate(exdAP=LifeTable(x=Age,mx=dAP,ax=0.5)$lt$ex,
         exdAC=LifeTable(x=Age,mx=dAC,ax=0.5)$lt$ex,
         exdA=LifeTable(x=Age,mx=dA,ax=0.5)$lt$ex,
         exdAPC=LifeTable(x=Age,mx=dAPC,ax=0.5)$lt$ex) %>%
  ungroup() %>%
  filter(Age==65) 

#read standard mortality rates
oldex <- read.table(file=file.path(data.dir, "pc_data.csv"), sep=",", header=TRUE) %>%
  mutate(Year=year) %>%
  filter(Age %in% ages, year %in% years) %>%
  group_by(country,sex,year) %>%
  mutate(ex=LifeTable(x=Age,mx=mx,ax=0.5)$lt$ex) %>%
  ungroup() %>%
  filter(Age==65) %>%
  mutate(country = recode(country,
                          "CHE" = "Switzerland",
                          "DNK" = "Denmark",
                          "FIN" = "Finland",
                          "FRATNP" = "France",
                          "GBRTENW" = "England and Wales",
                          "ITA" = "Italy",
                          "NLD" = "Netherlands",
                          "NOR" = "Norway",
                          "SWE" = "Sweden")) 

oldexind <- read.table(file=file.path(data.dir, "pc_data_si.csv"),sep=",", header = TRUE) %>%
  mutate(Year=floor(year), Age=floor(Age))%>%
  group_by(country,sex,year,Age) %>%
  summarise(Dx=sum(Dxi),Exp=sum(Exp)) %>%
  ungroup() %>%
  filter(Age %in% ages, year %in% years) %>%
  group_by(country,sex,year) %>%
  mutate(mx=Dx/Exp,
         exi=LifeTable(x=Age,mx=mx,ax=0.5)$lt$ex) %>%
  ungroup() %>%
  filter(Age==65) %>%
  mutate(country = recode(country,
                          "CHE" = "Switzerland",
                          "DNK" = "Denmark",
                          "FIN" = "Finland",
                          "FRATNP" = "France",
                          "GBRTENW" = "England and Wales",
                          "ITA" = "Italy",
                          "NLD" = "Netherlands",
                          "NOR" = "Norway",
                          "SWE" = "Sweden")) 

################################ together
ex65 <- newexTotal  %>%
  mutate(year=Per,country=countrysel, sex=sexsel) %>%
  select(year:sex,exdAP:exdAPC) %>%
  left_join(newexIndiv %>%
              mutate(year=Per,country=countrysel, sex=sexsel,
                     exdAPi=exdAP,exdAi=exdA,exdAPCi=exdAPC)%>%
              select(year:sex,exdAPi:exdAPCi)) %>%
  left_join(oldex %>%
              select(year,country,sex,ex)) %>%
  left_join(oldexind %>%
              select(year,country,sex,exi)) %>%
  pivot_longer(exdAP:exi, values_to = "ex", names_to = "type") %>%
  mutate(Type=ifelse((type=="ex"|type=="exi"),"Standard",type),
         Rates=ifelse((type=="ex"| type=="exdA"| type=="exdAP"| type=="exdAC"| type=="exdAPC"),"Total","Baseline")) 

## raw (unsmoothed) improvement series - not used by the plots below.
ex65impr <- ex65%>%
  arrange(country,sex,Type,Rates,year) %>%                                   ## [FIX] order before lag
  group_by(country,sex,Type,Rates) %>%
  summarise(year=year,impr=100*(ex-dplyr::lag(ex))/dplyr::lag(ex), .groups="drop") %>%  ## [FIX] dplyr::lag
  filter(!is.na(impr))


######plot appendix - standard, dAPC
femex <-  ggplot(ex65 %>% 
                   filter(sex == "Female", (type=="exdAPC"|type=="ex")) %>%
                   mutate(Type=ifelse(Type=="exdAPC","dAPC","Observed")),
                 aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line(size=0.7, cex=0.8)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "bottom") +
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10), 2019),
                     limits = c(1980, 2019)) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Female")


maleex <-  ggplot(ex65 %>% 
                    filter(sex == "Male", (type=="exdAPC"|type=="ex")) %>%
                    mutate(Type=ifelse(Type=="exdAPC","dAPC","Observed")),
                  aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line(size=0.7,cex=0.8)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "bottom") +
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10), 2019),
                     limits = c(1980, 2019)) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Male")


combined_plot_ex <- (femex + labs(subtitle = "Female")) +
  (maleex + labs(subtitle = "Male")) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(combined_plot_ex,filename=file.path(plot.dir, "exAPC.pdf"),width=17,height=10)




######plot appendix - baseline, dAPC
femexi <-  ggplot(ex65 %>% 
                   filter(sex == "Female", (type=="exdAPCi"|type=="exi")) %>%
                   mutate(Type=ifelse(Type=="exdAPCi","dAPC","Standard")),
                 aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line(size=0.7,cex=0.8)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "bottom") +
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10), 2019),
                     limits = c(1980, 2019)) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Female")


maleexi <-  ggplot(ex65 %>% 
                    filter(sex == "Male", (type=="exdAPCi"|type=="exi")) %>%
                    mutate(Type=ifelse(Type=="exdAPCi","dAPC","Standard")),
                  aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line(size=0.7,cex=0.8)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "bottom") +
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10), 2019),
                     limits = c(1980, 2019)) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Male")


combined_plot_ex <- (femexi + labs(subtitle = "Female")) +
  (maleexi + labs(subtitle = "Male")) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(combined_plot_ex,filename=file.path(plot.dir, "exAPCi.pdf"),width=17,height=10)


###### plot main - standard, dA, dAP , smooth standard first
  ex65smooth  <- ex65 %>%
    group_by(country, sex, type) %>%
    group_modify(~{
      dat <- arrange(.x, year)
      # Fit exactly as in the geom_smooth above
      fit <- gam(ex ~ s(year, k = 12, bs = "ps", sp = 3), data = dat)
      dat$ex_pspline <- predict(fit, newdata = dat) 
      dat}) %>%
    ungroup() %>%
    mutate(ex=ex_pspline) %>%
    select(-ex_pspline)
  
  imprex <- ex65smooth %>%
    arrange(country,sex,Type,Rates,year) %>%                                 ## [FIX] order before lag
    group_by(country,sex,Type, Rates) %>%
    summarise(year=year,impr=100*(ex-dplyr::lag(ex))/dplyr::lag(ex), .groups="drop") %>%  ## [FIX] dplyr::lag
    filter(!is.na(impr),Type!="exdAPC", Type!="exdAPCi") %>%
    mutate(Type=ifelse(Type=="exdAi","dA",Type),
           Type=ifelse(Type=="exdA","dA",Type),
           Type=ifelse(Type=="exdAPi","dAP",Type),
           Type=ifelse(Type=="exdAP","dAP",Type),
           Type=ifelse(Type=="exdACi","dAC",Type),
           Type=ifelse(Type=="exdACi","dAC",Type)) %>%
    mutate(Type  = factor(Type,  levels = c("Standard", "dAP", "dAC","dA")),
      Rates = factor(Rates, levels = c("Total", "Baseline")))


eximpr1fem <-  ggplot(imprex %>% filter(sex=="Female"), aes(x = year, y = impr, color=Type, linetype=Rates)) +
  geom_line(size=0.7)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12))+
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10),2019), limits = c(1980, 2019))+
  coord_cartesian(ylim = c(-0.2, 2.1)) +
  scale_y_continuous(breaks = seq(-0.2, 2.0, by = 0.2))+
  labs(x = "Year", y = "Change (in %)", color = "Type") +
  guides(linetype = guide_legend(override.aes = list(color = "black")))+
  theme(legend.position = "bottom")+
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.5)+ 
  labs(subtitle = "Female")

eximpr1male <- ggplot(imprex %>% filter(sex=="Male"), aes(x = year, y = impr, color=Type, linetype=Rates)) +
  geom_smooth(method = "gam", formula = y ~ s(x,k=10,bs="ps", sp=3), , se=FALSE)+
  facet_wrap(~ country) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12))+
  scale_x_continuous(breaks = c(seq(1980, 2010, by = 10),2019), limits = c(1980, 2019))+
  coord_cartesian(ylim = c(-0.2, 2.1)) +
  scale_y_continuous(breaks = seq(-0.2, 2.0, by = 0.2))+
  labs(x = "Year", y = "Change (in %)", color = "Type") +
  guides(linetype = guide_legend(override.aes = list(color = "black")))+
  theme(legend.position = "bottom")+
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.5)+ 
  labs(subtitle = "Male")



p1 <- eximpr1fem +
  labs(subtitle = "Female", color = "Type", linetype = "Rates") 
p2 <- eximpr1male +
  labs(subtitle = "Male", color = "Type", linetype = "Rates") 

combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")


ggsave(combined_plot,filename=file.path(plot.dir, "eximpr10w05s.pdf"),width=17,height=10)


###############################################################################
## compare cohort and period effects for total models, 1960-2019 vs 1980-2019
###############################################################################
## ---- read the extended-window (1960-2019) Total-rates effects --------------

coh1960   <- read.table(file = "cohorteffectbs1960.csv", sep = ",", header = TRUE) %>%
  distinct() %>% recode_country()
per1960   <- read.table(file = "pereffectbs1960.csv",    sep = ",", header = TRUE) %>%
  distinct() %>% recode_country()
drift1960 <- read.table(file = "driftbs1960.csv",        sep = ",", header = TRUE) %>%
  distinct() %>% recode_country()
drift1960 <- drift1960 %>% rename(drift = `exp.Est..`)

## ---- cohort effect: overlay the two windows -------------------------------
## cohTotal is the 1980-2019 Total cohort effect already read above.
coh_windows <- bind_rows(
  cohTotal %>% mutate(Window = "1980-2019"),
  coh1960  %>% mutate(Window = "1960-2019")
) %>%
  mutate(Window = factor(Window, levels = c("1980-2019", "1960-2019")))

## widen the x-axis so the extra (pre-1886) cohorts from the 1960 run are shown
coh_xmin   <- floor(min(coh_windows$Coh, na.rm = TRUE))
coh_breaks <- seq(1870, 1950, by = 20)

cohplot_windows <- ggplot(coh_windows,
                          aes(x = Coh, y = C.RR, color = sexsel, linetype = Window)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ countrysel) +
  scale_x_continuous(
    breaks = coh_breaks,
    limits = c(coh_xmin, 1954),
    labels = function(x) paste0(as.integer(x),
              "<br><span style='color:grey50;'>", as.integer(x) + 65, "</span>")) +
  scale_y_continuous(limits = c(0.5, 1.2), breaks = seq(0.5, 1.1, by = 0.1)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  labs(x = "Cohort / <span style='color:grey50;'>Year When 65</span>",
       y = "Relative Rate", color = "Sex", linetype = "Window") +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_markdown(lineheight = 0.95, margin = margin(t = 3)),
        axis.title.x = element_markdown(size = 12, margin = margin(t = 8)),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        panel.spacing.x = unit(1.2, "lines"))

ggsave(cohplot_windows, filename = "cohort_effect_1960_vs_1980.pdf", width = 10, height = 12)

###############################################################################
## Supplementary: drift-only restricted e-dagger and entropy (Hypothesis 1)
## Source: dA_edag_H.csv written by 5_LE_and_CI.R
###############################################################################

edagH <- read.table(file = "dA_edag_H.csv", sep = ",", header = TRUE) %>%
  distinct() %>%
  recode_country() %>%
  rename(country = countrysel, sex = sexsel)

edagH_panel <- function(yvar, ylab) {
  ggplot(edagH, aes(x = year, y = .data[[yvar]], color = sex, group = sex)) +
    geom_line(size = 0.7) +
    facet_wrap(~ country) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text   = element_text(size = 12, face = "bold"),
      axis.text    = element_text(size = 10),
      axis.title   = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "bottom") +
    scale_x_continuous(breaks = c(seq(1980, 2010, by = 10), 2019),
                       limits = c(1980, 2019)) +
    labs(x = "Year", y = ylab, color = "Sex", subtitle = NULL)
}

p1 <- edagH_panel("edag_dA", "Years") +
  labs(subtitle = "e-dagger")
p2 <- edagH_panel("H_dA", NULL) +
  labs(subtitle ="H")

combined_plot_edagH <- ggarrange(p1, p2, ncol = 2,
                                 common.legend = TRUE, legend = "bottom")

ggsave(combined_plot_edagH, filename = "edag_H_dA.pdf", width = 17, height = 10)
