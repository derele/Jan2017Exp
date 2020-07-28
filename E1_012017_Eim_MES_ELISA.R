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
colnames(E1)[1] <- "IFNy"
E1 <- dplyr::select(E1, IFNy)
setDT(E1, keep.rownames = TRUE)[]
colnames(E1)[1] <- "EH_ID"
colnames(E1_samples)[1] <- "EH_ID"
E1 <- merge(E1, E1_samples)
E1$OD <- NULL

E1 <- merge(E1, stab)
E1 <- merge(E1, M, by =c("EH_ID", "dpi.diss", "inf.strain") )
# E1 <- merge(E1, all.data, by = c("EH_ID", "dpi.diss", "inf.strain"))
# E1 <- E1[!duplicated(E1$IFNy),]

E1 <- dplyr::select(E1, EH_ID, IFNy, dpi.diss, inf.strain, dpi)
E1 <- dplyr::distinct(E1)
colnames(E1)[5] <- "dpi_count"
E1 <- merge(E1, all.data)

cor.test(E1$IFNy, E1$PH.delta, method = "kendall" )

ggplot(E1, aes( x = PH.delta, y = IFNy, color = inf.strain)) +   
  geom_point() +
  facet_wrap(~inf.strain)
as.numeric(E1$dpi.diss)
ggplot(E1, aes( x = dpi.diss, y = IFNy, color = inf.strain)) +   
  geom_point() +
  facet_wrap(~inf.strain)

write.csv(E1, "./Jan2017Exp/E1_012017_Eim_MES_ELISA1.csv")

