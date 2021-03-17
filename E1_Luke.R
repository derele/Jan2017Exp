library(httr)
library(RCurl)
library(Rmisc)
library(tidyr)
library(lme4)
library(ggplot2)
library(lmerTest)
library(modelr)
library(sjPlot)
library(ggpubr)

#### ### get the data --------------------------------------------

## for control: general experimental setup -------------------
stab <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv")
stab <- subset(stab, !stab$inf.strain%in%"EI70", drop = TRUE)
stab$inf.strain <- as.factor(as.character(stab$inf.strain))
names(stab)[names(stab)%in%"mouseID"] <- "EH_ID"

### IFNy MES
IFN_MES <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA.csv")
IFN_MES$X <- NULL

### IFNy SPL
IFN_SPL <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_SPL_ELISA.csv")
IFN_SPL$X <- NULL

### merge stab and IFN
IFN <- merge(IFN_MES, IFN_SPL, all = T)
stab <- merge(stab, IFN, all = T)

### spleen weight
spleen <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Spleen_weight.csv")

names(spleen)[names(spleen)%in%"mouseID"] <- "EH_ID"
spleen <- merge(spleen, stab)
spleen$dpi <- as.numeric(gsub("dpi|dip", "", spleen$dpi.diss))

## weight --------------------------------------------------------
weight <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Weight_Expe_Jan_2017.csv")

## weight was only obtained for mice dissected at or after 7dpi
## correcting the names
names(weight)[names(weight)%in%"Mouse.ID"] <- "EH_ID"

## EI70 infections were not followed up as likely double infections
weight <- weight[!weight$inf.strain%in%"EI70", ] 

## remove non meaningful columns
weight <- weight[, !colnames(weight)%in%c("Mouse.number", "Date.of.Birth")]

#percentage of weight per day
p.weight <- apply(weight[, grepl("^Day.*g", colnames(weight))], 2,
                  function (x) {
                    (x/weight$Day.1_.g)*100
                  })
#change gram to percentage weight column labels
colnames(p.weight) <- gsub(".g", ".p", colnames(p.weight))

#combine gram and percentage tables
weight <- cbind(weight, p.weight[, 2:ncol(p.weight)])

#reshape to long data format
weight.long <- reshape(weight,
                       direction = "long",
                       idvar = "EH_ID", ids = EH_ID,
                       varying = list(grep(".p$", colnames(weight), value=TRUE)),
                       timevar="dpi_of_perc",
                       v.names = "perc_of_dpi1", 
                       times = grep(".p$", colnames(weight), value=TRUE))
#make dpi numeric, keep just numbers
weight.long$dpi_of_perc <- as.numeric(gsub("Day\\.?(\\d+)_\\.p", "\\1",
                                           weight.long$dpi_of_perc))
#remove all Day.g columns
weight.long <- weight.long[, !grepl("^Day.", names(weight.long)) ]

### oocysts ---------------------------------------------------------
oocysts <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Clean_oocyst_data.csv")


#make total oocysts per g out of NB counts
oocysts$Total.oocysts.g <- ((oocysts$Count..8.Neubauer.squares. / 8)*
                              10000 * 2) / oocysts$used.in.flotation..g.

## correcting the names
names(oocysts)[names(oocysts)%in%c("Sample.ID", "dpi")] <- c("EH_ID", "dpi_count")
oocysts <- merge(oocysts, stab, all.y=TRUE)

all.data <- merge(oocysts[, c("EH_ID", "dpi_count", "Total.oocysts.g",
                              "dpi.diss", "inf.strain")],
                  weight.long,
                  by.x=c("EH_ID", "dpi_count", "dpi.diss", "inf.strain"),
                  by.y=c("EH_ID", "dpi_of_perc", "dpi.diss", "inf.strain"),
                  all=TRUE)

all.data$dpi.diss <- gsub("dip$", "dpi", all.data$dpi.diss)

## rename to the abreviations in the paper
levels(all.data$inf.strain) <- c("EfalL", "EfalW", "EferW", "Uninf", "EI70")

## we can observe some pattern:
## we have any measurements only for mice killed at or after 7dpi
## we have oocyst counts only for mice kileed at or after 7dpi

## in other words: 
## weight is not reported (measured) for mice killed at 3 and 5dpi
## oocyst counts are not repored (measured) for mice killed at 3, 5 and 7 dpi


############# --------------- qPCR for DNA -------------------
#set levels for infection strains
levels(stab$inf.strain) <- c("EfalL", "EfalW", "EferW", "Uninf")
#load in data
Rtissue <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Eimeria-NMRI_Relative%20quantification_clean.csv")

## only the means
RtMeans <- Rtissue[seq(1, nrow(Rtissue), by=2),
                   c("Sample", "Cq.Mean", "Cq.Mean.1")]
#naming and upper case
names(RtMeans) <- c("EH_ID", "Mouse_gDNA", "Eimeria_mDNA")
RtMeans$EH_ID <- toupper(RtMeans$EH_ID)

## LM0065 was measured twice with the same outcome
RtMeans <- RtMeans[!duplicated(RtMeans$EH_ID),]
RtMeans <- merge(RtMeans, stab, all=TRUE)
RtMeans$dpi.diss <- gsub("dip$", "dpi", RtMeans$dpi.diss)
RtMeans$dpi_count<- as.numeric(gsub("dpi", "", RtMeans$dpi.diss))

all.data <- merge(RtMeans, all.data, all=TRUE)

## remove the mice with no data whatsoever
all.data <- all.data[!is.na(all.data$dpi_count), ]
#Eim - Mouse = negative values are Eim positive
all.data$PH.delta <- all.data$Mouse_gDNA - all.data$Eimeria_mDNA
# Uninfected signal streght
max.neg <- max(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"],na.rm=TRUE)
# strongest Eim signal
min.neg <- min(all.data$PH.delta, na.rm=TRUE)
#baseline
LLD <- mean(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"], na.rm=TRUE)+
  2*(sd(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"], na.rm=TRUE))

no.log.lld <- 2^LLD
## per 100 transcripts (Dreisatz ;-))
(no.log.lld*1/(1-no.log.lld))*100

maxPHd <- max(all.data$PH.delta, na.rm=TRUE)
no.log.max.PHd <- 2^maxPHd
(no.log.max.PHd*1/(1-no.log.max.PHd))*100

max(all.data[all.data$inf.strain%in%"EferW", "PH.delta"], na.rm=TRUE)


### Histology ----------------------------------------------

hist <- read.csv("https://raw.githubusercontent.com/derele/Jan2017Exp/master/Histo.csv")
#add means and sums columns
hist$mMLS <- rowMeans(hist[, grepl("^MLS", colnames(hist))])
hist$sMLS <- rowSums(hist[, grepl("^MLS", colnames(hist))])

names(hist)[names(hist)%in%"Sample.ID"] <- "EH_ID"
all.data <- merge(all.data, hist, all=TRUE)

#order factors
all.data$dpi.diss <- factor(all.data$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))

# test strain effect on IFNy CEWE and MES + drop levels to remove EI70


IFNy.MES <- dplyr::select(all.data, EH_ID, dpi.diss, inf.strain, IFNy_MES)
IFNy.MES <- dplyr::distinct(IFNy.MES)
IFNy.MES <- na.omit(IFNy.MES)
IFNy.MES <- droplevels(IFNy.MES, exclude = "EI70")
IFNy.MES$inf.strain = factor(IFNy.MES$inf.strain, levels(IFNy.MES$inf.strain)[c(4,1:3)])
modIFNyMES <- lmer(IFNy_MES~inf.strain + (1|dpi.diss), data = IFNy.MES)
summary(modIFNyMES)
rand(modIFNyMES)
VarCorr(modIFNyMES)

MESgrid <- IFNy.MES %>% data_grid(inf.strain, dpi.diss, IFNy_MES)
MESgrid <- MESgrid %>% 
  add_predictions(modIFNyMES) 

ggplot(fortify(modIFNyMES), aes(dpi.diss, IFNy_MES, color = inf.strain, group = inf.strain)) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun=mean, geom="line")
# only Uninfected shows any effect (no strain difference)

# test strain effect on IFNy CEWE and SPL
IFNy.SPL <- dplyr::select(all.data, EH_ID, dpi.diss, inf.strain, IFNy_SPL)
IFNy.SPL <- dplyr::distinct(IFNy.SPL)
IFNy.SPL <- na.omit(IFNy.SPL)
IFNy.SPL <- droplevels(IFNy.SPL, exclude = "EI70")
IFNy.SPL$inf.strain = factor(IFNy.SPL$inf.strain, levels(IFNy.SPL$inf.strain)[c(4,1:3)])
modIFNySPL <- lmer(IFNy_SPL~inf.strain + (1|dpi.diss), data = IFNy.SPL)
summary(modIFNySPL)
rand(modIFNySPL)
VarCorr(modIFNySPL)
SPLgrid <- IFNy.SPL %>% data_grid(inf.strain, dpi.diss, IFNy_SPL)
SPLgrid <- SPLgrid %>%
  add_predictions(modIFNySPL)

ggplot(fortify(modIFNySPL), aes(dpi.diss, IFNy_SPL, color = inf.strain, group = inf.strain)) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun=mean, geom="line")

tab_model(modIFNyMES,modIFNySPL, 
          file="IFNtable_VS_Eflab(itercept).html",
          dv.labels=c("MES", "SPL"))

ggscatter(subset(all.data, !is.na(all.data$PH.delta)), x = "PH.delta", y = "IFNy_MES", add = "reg.line", color = "inf.strain") +
  facet_grid(~inf.strain) +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 510) +
  ggtitle("infection intensity effect on IFN-y abundance")

ggscatter(subset(all.data, !is.na(all.data$PH.delta)), x = "PH.delta", y = "IFNy_SPL", add = "reg.line", color = "inf.strain") +
  facet_grid(~inf.strain) +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 510) +
  ggtitle("infection intensity effect on IFN-y abundance")

# no differences

# sum of oocysts from dpi 4-9 vs IFN
oocysts <- dplyr::select(all.data, EH_ID, Total.oocysts.g, dpi_count, dpi.diss, IFNy_MES, IFNy_SPL)

# graph for presentation
ggplot(data = IFNy.MES, aes(x = dpi.diss, y = IFNy_MES, group = inf.strain, color = inf.strain)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  facet_grid(~inf.strain)

ggplot(data = IFNy.SPL, aes(x = dpi.diss, y = IFNy_SPL, group = inf.strain, color = inf.strain)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  facet_grid(~inf.strain)

ggplot(data = all.data, aes(x = PH.delta, y = IFNy_MES, color = inf.strain, group = inf.strain)) +
  geom_point() +
  facet_grid(~inf.strain)

############ - Gene expression data (spleen) -----------------------
#load data from raw GitHub
CXCL9.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_CXCL9.csv"
CXCL9.S <- read.csv(text = getURL(CXCL9.Surl), sep = ",")

IL10.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_IL10.csv"
IL10.S <- read.csv(text = getURL(IL10.Surl), sep = ",")

IL12.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_IL12.csv"
IL12.S <- read.csv(text = getURL(IL12.Surl), sep = ",")

IL6.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_IL6.csv"
IL6.S <- read.csv(text = getURL(IL6.Surl), sep = ",")

IFNg.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_INFg.csv"
IFNg.S <- read.csv(text = getURL(IFNg.Surl), sep = ",")

STAT6.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_STAT6.csv"
STAT6.S <- read.csv(text = getURL(STAT6.Surl), sep = ",")

TGFb.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_TGFb.csv"
TGFb.S <- read.csv(text = getURL(TGFb.Surl), sep = ",")

TNFa.Surl <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/GE_TNFa.csv"
TNFa.S <- read.csv(text = getURL(TNFa.Surl), sep = ",")

GeMeans.l <- list(CXCL9.S, IL10.S, IL12.S, IL6.S, IFNg.S, STAT6.S, TGFb.S, TNFa.S)

GeMeans.l <- lapply(GeMeans.l, function (data) {
  out <- data[!is.na(data$NE), c("Sample", "Gene", "NE")]
  out[!is.na(out$Gene), ]
})

GeMeans <- Reduce(rbind, GeMeans.l)

GeMeans$Sample <- toupper(GeMeans$Sample)
### Correction
## removing an empty row
GeMeans <- GeMeans[!GeMeans$Gene%in%"",]
GeMeans$Gene <- toupper(GeMeans$Gene)

## standard naming
names(GeMeans)[names(GeMeans)%in%"Sample"] <- "EH_ID"
## check uniqueness for genes / samples
nrow(unique(GeMeans)) ==  nrow(GeMeans)
nrow(unique(GeMeans[, c("EH_ID", "Gene")])) ==  nrow(GeMeans)
## Okay 456

## wide dateset for merging in overall table
GeMeans.wide <- reshape(GeMeans, timevar = "Gene", idvar = "EH_ID", direction = "wide")

M <- merge(GeMeans, stab, all=TRUE)
M$dpi <- as.numeric(gsub("dpi|dip", "", M$dpi.diss))
M$Gene <- as.character(M$Gene)
M$Gene[M$Gene == "IL6"] <- "IL-6"
M$Gene[M$Gene == "IL10"] <- "IL-10"
M$Gene[M$Gene == "IL12"] <- "IL-12"
M$Gene[M$Gene == "IFNG"] <- "IFN-G"
M$Gene[M$Gene == "INFG"] <- "IFN-G"
M$Gene[M$Gene == "TGFB"] <- "TGF-B"

M.wide <- merge(GeMeans.wide, stab, all=TRUE)

# add genes as factors to order
M$Gene_f = factor(M$Gene, levels=c('CXCL9','IL-6','IL-10','IL-12', "IFN-G", "TGF-B", "STAT6"))

# plot Spleen
CytokinesSP <- ggplot(subset(M, nchar(M$Gene)>2), aes(dpi, NE, color=inf.strain)) +
  geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  scale_x_continuous(breaks=c(3, 5, 7, 9, 11),
                     labels=c("3dpi", "5dpi", "7dip", "9dpi", "11dpi")) +
  facet_wrap(~Gene_f, scales="free_y", nrow=2)+
  scale_colour_brewer("infection\nisolate", palette = "Dark2") +
  scale_y_continuous("normalized mRNA expression")+
  theme_bw()

ggscatter(data = M, x = "NE", y = "IFNy_MES", add = "reg.line", color = "inf.strain") +
  facet_grid(Gene_f~inf.strain, scales = "free") +
  labs(y = "IFN-y (pg/mL)", x = "deltaCT = Mouse - Eimeria", color = "infection status", fill = "infection status") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  ylim(0, 510) +
  ggtitle("infection intensity effect on IFN-y abundance")

## Contrasting against Eflab
modCXCL9 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"CXCL9"))
summary(modCXCL9)

modIL10 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-10"))
summary(modIL10)

modIL12 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-12"))
summary(modIL12)

modIL6 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-6"))
summary(modIL6)

modIFNG <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IFN-G"))
summary(modIFNG)

modSTAT6 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"STAT6"))
summary(modSTAT6)

modTGFB <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TGF-B"))
summary(modTGFB)

modTNFA <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TNFA"))
summary(modTNFA)

tab_model(modCXCL9, modIL10, modIL12, modIL6,
          modIFNG, modSTAT6, modTGFB, 
          file="SPtable_VS_Eflab(itercept).html",
          dv.labels=c("CXCL9", "IL10", "IL12", "IL6",
                      "IFNG", "STAT6", "TGFB"))

# contrast against uninfected
M$inf.strain = factor(M$inf.strain, levels(M$inf.strain)[c(4,1:3)])

modCXCL9.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"CXCL9"))
summary(modCXCL9.cu)

modIL10.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-10"))
summary(modIL10.cu)

modIL12.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-12"))
summary(modIL12.cu)

modIL6.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL-6"))
summary(modIL6.cu)

modIFNG.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IFN-G"))
summary(modIFNG.cu)

modSTAT6.cu <- lmer(NE~inf.strain  +(1|dpi.diss), data=subset(M, M$Gene%in%"STAT6"))
summary(modSTAT6.cu)

modTGFB.cu <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TGF-B"))
summary(modTGFB.cu)

ggplot(data = all.data, aes(x = PH.delta, y = IFNy_SPL, color = inf.strain, group = inf.strain)) +
  geom_point() +
  facet_grid(~inf.strain)
