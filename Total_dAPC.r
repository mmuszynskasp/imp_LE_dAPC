#################################################################################################################################
####Aim: fit the dAPC model, main model with 10 knots and weighted trend

rm(list=ls())

library(tidyr)
library(dplyr)
library(MortalityLaws)
library(Epi)
library(HMDHFDplus)
library(mgcv)
library(splines)
library(patchwork)

models.dir <- #your directory to save the model results data
out.dir <- #your directory with data prepared in the last step


setwd(out.dir)
lexis <- read.table(file="lexis.csv", sep=",", header=TRUE)

countries <- unique(lexis$country)
sexes <- unique(lexis$sex)


########the first model to set up the files and later append the results
i=1
j=1

countrysel <- countries[i]
sexsel <- sexes[j]

mylexis <- lexis %>%
  mutate(D=round(Dx), Y=round(Exp),
         A = Age, P = Year) %>%
  filter(country==countrysel, sex==sexsel)

write.table(mylexis,file="mydataapc.csv",sep=",")
mydataapc <- as.data.frame(read.table(file="mydataapc.csv",sep=",",header=TRUE))

mymodel2 <- apc.fit(data=mydataapc,npar=c(A=10,P=10,C=10), model="bs",  parm="APC", dr.extr="Y",ref.p=1980+1/3)

setwd(models.dir)
cohorteff <- cbind(mymodel2$Coh[,1:2],countrysel,sexsel)
write.table(cohorteff, file="cohorteffectbs10w.csv",sep=",",row.names = FALSE)

pereff <- cbind(mymodel2$Per[,1:2],countrysel,sexsel)
write.table(pereff, file="pereffectbs10w.csv",sep=",",row.names = FALSE)


ageeff <- cbind(mymodel2$Age[,1:2],countrysel,sexsel)
write.table(ageeff, file="ageeffectbs10w.csv",sep=",",row.names = FALSE)

drift <- mymodel2$Drift[1,1]
mydrift <- cbind(drift,countrysel,sexsel)
write.table(mydrift, file="driftbs10w.csv",sep=",",row.names = FALSE)

for (i in 1:length(countries)){  #make the full loop again and then remove the duplicates when reading the data
  for (j in 1:length(sexes)){
    countrysel <- countries[i]
    sexsel <- sexes[j]
    
    mylexis <- lexis %>%
      mutate(D=round(Dx), Y=round(Exp),
             A = Age, P = Year) %>%
      filter(country==countrysel, sex==sexsel)
    
    write.table(mylexis,file="mydataapc.csv",sep=",")
    mydataapc <- as.data.frame(read.table(file="mydataapc.csv",sep=",",header=TRUE))
    
    mymodel2 <- apc.fit(data=mydataapc,npar=c(A=10,P=10,C=10), model="bs",  parm="APC",dr.extr="Y",ref.p=1980+1/3)

    setwd(models.dir)
    cohorteff <- cbind(mymodel2$Coh[,1:2],countrysel,sexsel)
    write.table(cohorteff, file="cohorteffectbs10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    pereff <- cbind(mymodel2$Per[,1:2],countrysel,sexsel)
    write.table(pereff, file="pereffectbs10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    
    ageeff <- cbind(mymodel2$Age[,1:2],countrysel,sexsel)
    write.table(ageeff, file="ageeffectbs10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
    drift <- mymodel2$Drift[1,1]
    mydrift <- cbind(drift,countrysel,sexsel)
    write.table(mydrift, file="driftbs10w.csv",sep=",",row.names = FALSE, col.names=FALSE, append=TRUE)
    
  }    
} 

