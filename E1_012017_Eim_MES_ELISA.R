# E1 012017 Eim ELISA cleanup and processing into clean table
library(httr)
library(RCurl)
library(Rmisc)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

E1_std <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA1_std.csv"))
E1_samples <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA1_samples.csv"))

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E1_std)
plot(model1)

E1<-ED(model1, E1_samples$OD, type="absolute", display=F)
row.names(E1) <- E1_samples$Sample

points(y=E1_samples$OD,x=E1[,1],col="lightblue",pch=19,cex=2)
text(y =E1_samples$OD, x = E1[,1], labels=E1_samples$Sample, data=E1, cex=0.9, font=2)

E1 <- data.frame(E1)
colnames(E1)[1] <- "IFNy_CEWE"
E1 <- dplyr::select(E1, IFNy_CEWE)
setDT(E1, keep.rownames = TRUE)[]
colnames(E1)[1] <- "EH_ID"
colnames(E1_samples)[1] <- "EH_ID"
E1 <- merge(E1, E1_samples)
E1$OD <- NULL

write.csv(E1, "./Jan2017Exp/E1_012017_Eim_MES_ELISA1.csv")

E2_std <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA2_std.csv"))
E2_samples <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA2_samples.csv"))

model2<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E2_std)
plot(model2)

E2<-ED(model2, E2_samples$OD, type="absolute", display=F)
row.names(E2) <- E2_samples$Sample

points(y=E2_samples$OD,x=E2[,1],col="lightblue",pch=19,cex=2)
text(y =E2_samples$OD, x = E2[,1], labels=E2_samples$Sample, data=E2, cex=0.9, font=2)

E2 <- data.frame(E2)
colnames(E2)[1] <- "IFNy_CEWE"
E2 <- dplyr::select(E2, IFNy_CEWE)
setDT(E2, keep.rownames = TRUE)[]
colnames(E2)[1] <- "EH_ID"
colnames(E2_samples)[1] <- "EH_ID"
E2 <- merge(E2, E2_samples)
E2$OD <- NULL

E1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA1.csv"))
E1$X <- NULL
E <- rbind(E1, E2)

write.csv(E, "./Jan2017Exp/E1_012017_Eim_MES_ELISA.csv")

