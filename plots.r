#################################################################################################################################
####Aim: present the results of the dAPC models
rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)
library(stringr)
library(purrr)
library(zoo)
library(tibble)
library(mgcv)
library(patchwork)
library(MortalityLaws)

models.dir <- "your directory"
out.dir <- "your directory"
plot.dir <- "your directory"


ourcohorts <- 1885:1954
years <- 1980:2019 
ages <- 65:94

#########read in data from both models
setwd(models.dir)

cohTotal <- read.table(file="cohorteffectbs10w.csv",sep=",",header=TRUE) %>%
  distinct()
cohIndiv <- read.table(file="cohorteffectindiv10w.csv",sep=",",header=TRUE)%>%
  distinct()
perTotal <- read.table(file="pereffectbs10w.csv",sep=",",header=TRUE)%>%
  distinct()
perIndiv <- read.table(file="pereffectindiv10w.csv",sep=",",header=TRUE)%>%
  distinct()
ageTotal <- read.table(file="ageeffectbs10w.csv",sep=",",header=TRUE)%>%
  distinct()
ageIndiv <- read.table(file="ageeffectindiv10w.csv",sep=",",header=TRUE)%>%
  distinct()
driftTotal <- read.table(file="driftbs10w.csv",sep=",", header=TRUE) %>%
  distinct()
driftIndiv <- read.table(file="driftindiv10w.csv",sep=",", header=TRUE)%>%
  distinct()

anovaTotal <- read.table(file="myanovabs10w.csv",sep=",",header=TRUE) %>% 
  rename("countrysel"="country") #for replacement of country names in the next step
anovaIndiv <- read.table(file="myanovaindiv10w.csv",sep=",",header=TRUE)%>% 
  rename("countrysel"="country")#for replacement of country names in the next step

###recode country names for the plots
cntry_map <- c(CHE = "Switzerland", DNK = "Denmark",FIN = "Finland",FRATNP = "France",GBRTENW = "England and Wales",  
               ITA = "Italy",NLD = "Netherlands", NOR = "Norway",SWE = "Sweden")

recode_country <- function(df) {
  df %>%
    mutate(countrysel = recode(countrysel, !!!as.list(cntry_map),
                               .default = countrysel))}

objs <- c("cohTotal","cohIndiv","perTotal","perIndiv","ageTotal","ageIndiv","driftTotal","driftIndiv","anovaTotal","anovaIndiv")

updated <- lapply(mget(objs, inherits = TRUE), recode_country)
list2env(updated, .GlobalEnv)

#######################plots
## Contribution of effects
totest <- anovaTotal %>%
  unique()%>%
  mutate(Type="Total") %>%
  add_row(anovaIndiv %>%
            unique()%>%
            mutate(Type="Baseline")) %>%
  pivot_wider(names_from = model, values_from = dev) %>%
  mutate(country=countrysel,
         reddrift= Age - `Age-drift`,
         redP = `Age-drift`-`Age-Period`,
         redC = `Age-Period` - `Age-Period-Cohort`,
         total= Age-`Age-Period-Cohort`,
         Drift=100*reddrift/total,
         Period=100*redP/total,
         Cohort=100*redC/total) %>%
  select(country,sex,Cohort,Period,Type) %>%
  pivot_longer(Cohort:Period, names_to = "Effect", values_to= "Contribution")


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


setwd(plot.dir)
ggsave(anovaplot,filename="anovaplot.pdf",width=11,height=14)

### table
tab_anova <- left_join(totest %>% filter(sex=="Female") %>%
  pivot_wider(names_from = c(Effect,Type), values_from = Contribution, names_glue = "{Effect}|{Type}") %>%
  arrange(country) %>%
  select(country,`Cohort|Total`,`Cohort|Baseline`,`Period|Total`,`Period|Baseline`) %>%
    mutate(`Drift|Total`=100-`Cohort|Total`-`Period|Total`,`Drift|Baseline`=100-`Cohort|Baseline`-`Period|Baseline`),
  totest %>% filter(sex=="Male") %>%
    pivot_wider(names_from = c(Effect,Type), values_from = Contribution, names_glue = "{Effect}|{Type}") %>%
    select(country,`Cohort|Total`,`Cohort|Baseline`,`Period|Total`,`Period|Baseline`)%>%
    mutate(`Drift|Total`=100-`Cohort|Total`-`Period|Total`,`Drift|Baseline`=100-`Cohort|Baseline`-`Period|Baseline`), by="country")

print(xtable(tab_anova, digits=1), include.rownames = FALSE)
tab_anova_aver <- colMeans(tab_anova[,-1])


###############################################
#cohort effect
#we take cohort with maximum effect and cohort with minium (some limits to min for nice graphs)
coheffbsmin <- cohTotal %>%
  filter(Coh>1910, Coh<1950) %>%
  group_by(countrysel,sexsel) %>%
  dplyr::summarise(min=min(log(C.RR))) 

coheffbsmax <- cohTotal %>%
  filter(Coh>1910, Coh<1945) %>%  # 1945 is the correction for the last observations in Finland
  group_by(countrysel,sexsel) %>%
  dplyr::summarise(max=max(log(C.RR)))

minmaxyear <- coheffbsmin %>%
  left_join(cohTotal) %>%
  filter(min==log(C.RR)) %>%
  rename("Cohmin"="Coh") %>%
  mutate(Cohmin=ifelse(Cohmin==1949,NA,Cohmin)) %>%
  left_join(coheffbsmax %>%
              left_join(cohTotal) %>%
              filter(max==log(C.RR)) %>%
              rename("Cohmax"="Coh") %>%
              select(-C.RR)) %>%
  mutate(yearmax=Cohmax+65, yearmin=Cohmin+65)


cohplot <- ggplot(cohTotal, aes(x = Coh, y = C.RR, color = sexsel)) +
  geom_line(size = 1) +
  geom_vline(data = minmaxyear, aes(xintercept = Cohmin, color = sexsel), linetype = "dashed", size = 0.5) +
  geom_vline(data = minmaxyear, aes(xintercept = Cohmax, color = sexsel), linetype = "dashed", size = 0.5) +
  # Add year labels at the bottom with different positioning based on sexsel
  geom_text(data = minmaxyear, aes(x = Cohmin, y = 0.55, label = Cohmin, color = sexsel), 
            size = 3, 
            vjust = ifelse(minmaxyear$sexsel == "Male", -1, 1),  # Top for Male, Bottom for Female
            hjust = 1.2) +  # Center the labels horizontally
  geom_text(data = minmaxyear, aes(x = Cohmax, y = 0.55, label = Cohmax, color = sexsel), 
            size = 3, 
            vjust = ifelse(minmaxyear$sexsel == "Male", -1, 1),  # Top for Male, Bottom for Female
            hjust = 1.2) +
  facet_wrap(~ countrysel) +
  scale_x_continuous(breaks = c(1886, seq(1900, 1940, by = 10), 1954), limits = c(1886, 1954)) +
  scale_y_continuous(limits = c(0.5, 1.15), breaks = seq(0.5, 1.1, by = 0.1)) + 
  labs(x = "Cohort", y = "Relative Rate", color = "Sex") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12))+                
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")

setwd(plot.dir)
ggsave(cohplot,filename="cohlpt10w.pdf",width=10,height=10)


###### check when RR=1 for the first time - only in text
first_over1 <- cohTotal %>%
  group_by(countrysel, sexsel) %>%
  arrange(Coh, .by_group = TRUE) %>%
  summarise(first_cohort = {over <- .data[["C.RR"]] > 1
      if (any(over, na.rm = TRUE)) min(.data[["Coh"]][over], na.rm = TRUE) else NA_real_},.groups = "drop") 

####Baseline
coheffbsminI <- cohIndiv %>%
  filter(Coh>1910, Coh<1950) %>%
  group_by(countrysel,sexsel) %>%
  dplyr::summarise(min=min(log(C.RR))) 

coheffbsmaxI <- cohIndiv %>%
  filter(Coh>1910, Coh<1945) %>%  
  group_by(countrysel,sexsel) %>%
  dplyr::summarise(max=max(log(C.RR)))

minmaxyearI <- coheffbsminI %>%
  left_join(cohIndiv) %>%
  filter(min==log(C.RR)) %>%
  rename("Cohmin"="Coh") %>%
  mutate(Cohmin=ifelse(Cohmin==1949,NA,Cohmin)) %>%
  left_join(coheffbsmaxI %>%
              left_join(cohIndiv) %>%
              filter(max==log(C.RR)) %>%
              rename("Cohmax"="Coh") %>%
              select(-C.RR)) %>%
  mutate(yearmax=Cohmax+65, yearmin=Cohmin+65)


cohplotI <- ggplot(cohIndiv, aes(x = Coh, y = C.RR, color = sexsel)) +
  geom_line(size = 1) +
  geom_vline(data = minmaxyearI, aes(xintercept = Cohmin, color = sexsel), linetype = "dashed", size = 0.5) +
  geom_vline(data = minmaxyearI, aes(xintercept = Cohmax, color = sexsel), linetype = "dashed", size = 0.5) +
  # Add year labels at the bottom with different positioning based on sexsel
  geom_text(data = minmaxyearI, aes(x = Cohmin, y = 0.55, label = Cohmin, color = sexsel), 
            size = 3, 
            vjust = ifelse(minmaxyearI$sexsel == "Male", -1, 1),  # Top for Male, Bottom for Female
            hjust = 1.2) +  # Center the labels horizontally
  geom_text(data = minmaxyearI, aes(x = Cohmax, y = 0.55, label = Cohmax, color = sexsel), 
            size = 3, 
            vjust = ifelse(minmaxyearI$sexsel == "Male", -1, 1),  # Top for Male, Bottom for Female
            hjust = 1.2) +
  facet_wrap(~ countrysel) +
  scale_x_continuous(breaks = c(1886, seq(1900, 1940, by = 10), 1954), limits = c(1886, 1954)) +
  scale_y_continuous(limits = c(0.5, 1.15), breaks = seq(0.5, 1.1, by = 0.1)) + 
  labs(x = "Cohort", y = "Relative Rate", color = "Sex") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12))+                
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")

setwd(plot.dir)
ggsave(cohplotI,filename="cohlptindiv10w.pdf",width=10,height=10)


#####################################################################################################################################################
##### period effect combined with the period effect of the standard APC model
pereffbs <- perTotal %>%
    left_join(driftTotal) %>%
  group_by(countrysel,sexsel)%>%
  mutate(drifty=row_number()-1,
         driftobs=(drift^0.5)^drifty,
         newP.RR=P.RR/driftobs) %>%
  mutate(Per=floor(Per)) %>%
  group_by(Per,countrysel,sexsel) %>%       ####necessary because Lexis triangles
  summarise(P.RR=mean(as.numeric(P.RR)), 
            newP.RR=mean(as.numeric(newP.RR))) %>%
  ungroup() %>%  
  mutate(Rates="Total") %>%
  add_row(perIndiv %>%
            left_join(driftIndiv) %>%
            group_by(countrysel,sexsel)%>%
            mutate(drifty=row_number()-1,
                   driftobs=(drift^0.5)^drifty,
                   newP.RR=P.RR/driftobs) %>%
            mutate(Per=floor(Per)) %>%
            group_by(Per,countrysel,sexsel) %>%       ####necessary because Lexis triangles
            summarise(P.RR=mean(as.numeric(P.RR)), 
                      newP.RR=mean(as.numeric(newP.RR))) %>%
            ungroup() %>%
            mutate(Rates="Baseline")) %>%
  group_by(Rates,countrysel,sexsel) %>%
  mutate(change= -100*(P.RR-lag(P.RR))/lag(P.RR)) %>%
  mutate(Rates = factor(Rates, levels = c("Total", "Baseline"))) 



####plot without the drift, annual RR
perplotnew <- ggplot(pereffbs, aes(x = Per, y = newP.RR, color = sexsel, linetype=Rates)) +
  geom_line(size=1)+
  #geom_smooth(method = "gam", formula = y ~ s(x,k=10,bs="ps", sp=3), size=1, se=FALSE)+
  scale_color_manual(values = c("black", "red"))+
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
  #  geom_vline(xintercept = 2010, linetype = "dashed", color = "black", size = 0.5)+
  scale_x_continuous(breaks = c(seq(1980, 2015, by = 10), 2019), limits = c(1980, 2019))+
 # scale_y_continuous(breaks = c(seq(0, 1, by = 0.1)), limits = c(0, 1))+
  # scale_y_continuous(breaks = c(seq(-4.5, 0.5, by = 1),0), limits = c(-4.5, 0.5))
  theme(legend.position = "bottom")+ 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")

setwd(plot.dir)
ggsave(perplotnew,filename="perplotnew.pdf",width=10,height=10)


###################table for drift
drift <- driftTotal%>%
  mutate(drift=drift^0.5) %>%
  mutate(Rates="Total") %>%
  add_row(driftIndiv %>%
          mutate(Rates="Baseline"))

tabledrift <- drift %>%
  filter(sexsel=="Female") %>%
  select(-sexsel) %>%
  pivot_wider(names_from = Rates, values_from = drift) %>%
  left_join(drift %>%
              filter(sexsel=="Male") %>%
              select(-sexsel) %>%
              mutate(Rates=ifelse(Rates=="Total","Totalm","Baselinem")) %>%
              pivot_wider(names_from = Rates, values_from = drift)) %>%
  arrange(countrysel)

library(xtable)
print(xtable(tabledrift, digits=3), include.rownames = FALSE)

averdrift <- colMeans(tabledrift[,-1])
print(xtable(t(averdrift), digits=3), include.rownames = FALSE)

#####age effect, both types
ageeffbs <- ageTotal  %>%
  distinct() %>%
  mutate(Type="Total") %>%
  add_row(ageIndiv  %>%
            distinct() %>%
            mutate(Type="Baseline"))



ageplot <- ggplot(ageeffbs, aes(x = Age, y = log(Rate), color = sexsel, linetype=Type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("black", "red"))+
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
#scale_y_continuous(breaks = c(seq(-5, 0, by = 1)), limits = c(-5, 0))


setwd(plot.dir)
ggsave(ageplot,filename="ageplot.pdf",width=10,height=10)

#####################################################################################################
############################ LE plots

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
         dA=Rate*(drift^(Per-1980)),
         dAPC=Rate*C.RR*P.RR) %>%
  group_by(countrysel,sexsel,Per) %>%
  mutate(exdAP=LifeTable(x=Age,mx=dAP,ax=0.5)$lt$ex,
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
         dA=Rate*(drift^(Per-1980)),
         dAPC=Rate*C.RR*P.RR) %>%
  group_by(countrysel,sexsel,Per) %>%
  mutate(exdAP=LifeTable(x=Age,mx=dAP,ax=0.5)$lt$ex,
         exdA=LifeTable(x=Age,mx=dA,ax=0.5)$lt$ex,
         exdAPC=LifeTable(x=Age,mx=dAPC,ax=0.5)$lt$ex) %>%
  ungroup() %>%
  filter(Age==65) 

#read standard mortality rates
setwd(out.dir)
oldex <- read.table(file="pc_data.csv", sep=",", header=TRUE) %>%
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

oldexind <- read.table(file="pc_data_si.csv",sep=",", header = TRUE) %>%
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
         Rates=ifelse((type=="ex"| type=="exdA"| type=="exdAP"| type=="exdAPC"),"Total","Baseline")) 
 
ex65impr <- ex65%>%
  group_by(country,sex,Type,Rates) %>%
  summarise(year=year,impr=100*(ex-lag(ex))/lag(ex)) %>%
  filter(!is.na(impr))

#mutate(Rates = factor(Rates, levels = c("Total", "Baseline")),
#       Type = factor(Type, levels = c("Standard", "dAP")))


######plot appendix - standard, dAPC
femex <-  ggplot(ex65 %>% 
                   filter(sex == "Female", (type=="exdAPC"|type=="ex")) %>%
                   mutate(Type=ifelse(Type=="exdAPC","dAPC","Standard")),
                 aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line()+
  # geom_smooth(method = "gam",formula = y ~ s(x, k = 10, bs = "ps", sp = 3),se = FALSE,linewidth = 1) +
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
  # coord_cartesian(ylim = c(-0.2, 2)) +
  #  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Female")


maleex <-  ggplot(ex65 %>% 
                    filter(sex == "Male", (type=="exdAPC"|type=="ex")) %>%
                    mutate(Type=ifelse(Type=="exdAPC","dAPC","Standard")),
                  aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line()+
  # geom_smooth(method = "gam",formula = y ~ s(x, k = 10, bs = "ps", sp = 3),se = FALSE,linewidth = 1) +
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
  # coord_cartesian(ylim = c(-0.2, 2)) +
  #  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Male")


combined_plot_ex <- (femex + labs(subtitle = "Female")) +
  (maleex + labs(subtitle = "Male")) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

setwd(plot.dir)
ggsave(combined_plot_ex,filename="exAPC.pdf",width=17,height=10)




######plot appendix - baseline, dAPC
femexi <-  ggplot(ex65 %>% 
                   filter(sex == "Female", (type=="exdAPCi"|type=="exi")) %>%
                   mutate(Type=ifelse(Type=="exdAPCi","dAPC","Standard")),
                 aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line()+
  # geom_smooth(method = "gam",formula = y ~ s(x, k = 10, bs = "ps", sp = 3),se = FALSE,linewidth = 1) +
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
  # coord_cartesian(ylim = c(-0.2, 2)) +
  #  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Female")


maleexi <-  ggplot(ex65 %>% 
                    filter(sex == "Male", (type=="exdAPCi"|type=="exi")) %>%
                    mutate(Type=ifelse(Type=="exdAPCi","dAPC","Standard")),
                  aes(x = year, y = ex, color = Type, group = Type)) +
  geom_line()+
  # geom_smooth(method = "gam",formula = y ~ s(x, k = 10, bs = "ps", sp = 3),se = FALSE,linewidth = 1) +
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
  # coord_cartesian(ylim = c(-0.2, 2)) +
  #  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(x = "Year", y = "Years", color = "Type", subtitle = "Male")


combined_plot_ex <- (femexi + labs(subtitle = "Female")) +
  (maleexi + labs(subtitle = "Male")) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

setwd(plot.dir)
ggsave(combined_plot_ex,filename="exAPCi.pdf",width=17,height=10)


###### plot main - standard, dA, dAP , smooth standard first
  ex65smooth  <- ex65 %>%
    group_by(country, sex, type) %>%
    group_modify(~{
      dat <- arrange(.x, year)
      # Fit exactly as in your geom_smooth
      fit <- gam(ex ~ s(year, k = 12, bs = "ps", sp = 3), data = dat)
      dat$ex_pspline <- predict(fit, newdata = dat) 
      dat}) %>%
    ungroup() %>%
    mutate(ex=ex_pspline) %>%
    select(-ex_pspline)
  
  imprex <- ex65smooth %>%
    group_by(country,sex,Type, Rates) %>%
    summarise(year=year,impr=100*(ex-lag(ex))/lag(ex)) %>%
    filter(!is.na(impr),Type!="exdAPC", Type!="exdAPCi") %>%
    mutate(Type=ifelse(Type=="exdAi","dA",Type),
           Type=ifelse(Type=="exdA","dA",Type),
           Type=ifelse(Type=="exdAPi","dAP",Type),
           Type=ifelse(Type=="exdAP","dAP",Type)) %>%
    mutate(Type  = factor(Type,  levels = c("Standard", "dAP", "dA")),
      Rates = factor(Rates, levels = c("Total", "Baseline")))


eximpr1fem <-  ggplot(imprex %>% filter(sex=="Female"), aes(x = year, y = impr, color=Type, linetype=Rates)) +
  geom_line(size=1)+
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
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5)+ 
  labs(subtitle = "Female")

eximpr1male <- ggplot(imprex %>% filter(sex=="Male"), aes(x = year, y = impr, color=Type, linetype=Rates)) +
  geom_smooth(method = "gam", formula = y ~ s(x,k=10,bs="ps", sp=3), size=1, se=FALSE)+
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
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5)+ 
  labs(subtitle = "Male")



p1 <- eximpr1fem +
  labs(subtitle = "Female", color = "Type", linetype = "Rates") 
p2 <- eximpr1male +
  labs(subtitle = "Male", color = "Type", linetype = "Rates") 

combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")


setwd(plot.dir)
ggsave(combined_plot,filename="eximpr10w05s.pdf",width=17,height=10)



