#################################################################################################################################
####Aim: fit the dAPC model for baseline mortality rates
#####    model with 10 knots and weighted trend

rm(list=ls())

library(tidyr)
library(dplyr)
library(Epi)
library(mgcv)
library(splines)

models.dir <- "your directory"
out.dir <-"your directory"
plot.dir <- "your directory"


ourcohorts <- 1885:1954
years <- 1980:2019 
ages <- 65:94

setwd(out.dir)

pc_data_si <- read.table(file="lexis.csv", sep=",", header=TRUE) %>%
  mutate(Cohort2=Cohort, Cohort=floor(Cohort), year=floor(Year), age=floor(Age)) %>%
  left_join(read.table(file=paste(models.dir,"multip.csv",sep=""),sep=",", header = TRUE) %>%
              rename("year"="Year", "age"="Age") %>%
              mutate(sex=ifelse(sex=="fem","Female","Male"))) %>%
  mutate(Dxi=Dx*(multip^(-0.25)),  ##### adjustment for the frailty variance=0.25
         D=round(Dxi), Y=round(Exp),
         A = Age, P = Year)
  
setwd(out.dir)
write.table(pc_data_si, file="pc_data_si.csv",sep=",", row.names = TRUE)

individual <- read.table(file="pc_data_si.csv",sep=",", header = TRUE)

#################### fit the models

countries <- unique(individual$country)
sexes <- unique(individual$sex)

i=1
j=1
countrysel <- countries[i]
sexsel <- sexes[j]

myindividual <- individual  %>%
  filter(country==countrysel, sex==sexsel)


write.table(myindividual,file="mydataapc.csv",sep=",")
mydataapc <- as.data.frame(read.table(file="mydataapc.csv",sep=",",header=TRUE))

mymodel2 <- apc.fit(data=mydataapc,npar=c(A=10,P=10,C=10), model="bs",  parm="APC", dr.extr="Y",ref.p=1980+1/3)

setwd(models.dir)
cohorteff <- cbind(mymodel2$Coh[,1:2],countrysel,sexsel)
write.table(cohorteff, file="cohorteffectindiv10w.csv",sep=",",row.names = FALSE)

pereff <- cbind(mymodel2$Per[,1:2],countrysel,sexsel)
write.table(pereff, file="pereffectindiv10w.csv",sep=",",row.names = FALSE)


ageeff <- cbind(mymodel2$Age[,1:2],countrysel,sexsel)
write.table(ageeff, file="ageeffectindiv10w.csv",sep=",",row.names = FALSE)

drift <- mymodel2$Drift[1,1]
mydrift <- cbind(drift,countrysel,sexsel)
write.table(mydrift, file="driftindiv10w.csv",sep=",",row.names = FALSE)

myanova <- cbind(mymodel2$Anova$Model, mymodel2$Anova$`Mod. dev.`, countrysel,sexsel)
colnames(myanova) <- c("model","dev","country", "sex")
write.table(myanova, file="myanovaindiv10w.csv",sep=",",row.names = FALSE)

for (i in 1:length(countries)){
  for (j in 1:length(sexes)){
    countrysel <- countries[i]
    sexsel <- sexes[j]
    
    myindividual <- individual  %>%
      filter(country==countrysel, sex==sexsel)
    
    write.table(myindividual,file="mydataapc.csv",sep=",")
    mydataapc <- as.data.frame(read.table(file="mydataapc.csv",sep=",",header=TRUE))
    
    mymodel2 <- apc.fit(data=mydataapc,npar=c(A=10,P=10,C=10), model="bs",  parm="APC", dr.extr="Y",ref.p=1980+1/3)
    
    setwd(models.dir)
    cohorteff <- cbind(mymodel2$Coh[,1:2],countrysel,sexsel)
    write.table(cohorteff, file="cohorteffectindiv10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    pereff <- cbind(mymodel2$Per[,1:2],countrysel,sexsel)
    write.table(pereff, file="pereffectindiv10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    
    ageeff <- cbind(mymodel2$Age[,1:2],countrysel,sexsel)
    write.table(ageeff, file="ageeffectindiv10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    drift <- mymodel2$Drift[1,1]
    mydrift <- cbind(drift,countrysel,sexsel)
    write.table(mydrift, file="driftindiv10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    myanova <- cbind(mymodel2$Anova$Model, mymodel2$Anova$`Mod. dev.`, countrysel,sexsel)
    colnames(myanova) <- c("model","dev","country", "sex")
    write.table(myanova, file="myanovaindiv10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
  }
}


