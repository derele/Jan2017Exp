library(httr)
library(RCurl)
library(Rmisc)
library(tidyr)

#### ### get the data --------------------------------------------

## for control: general experimental setup -------------------
stabURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Experiment_Table_raw_NMRI_Jan2017.csv"
stab <- read.csv(text = getURL(stabURL))
stab <- subset(stab, !stab$inf.strain%in%"EI70", drop = TRUE)
stab$inf.strain <- as.factor(as.character(stab$inf.strain))
names(stab)[names(stab)%in%"mouseID"] <- "EH_ID"

### IFNy MES
IFN_MES <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_MES_ELISA.csv"))
IFN_MES$X <- NULL

### IFNy SPL
IFN_SPL <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Jan2017Exp/master/E1_012017_Eim_SPL_ELISA.csv"))
IFN_SPL$X <- NULL

### merge stab and IFN
IFN <- merge(IFN_MES, IFN_SPL, all = T)
stab <- merge(stab, IFN, all = T)

### spleen weight
spleenURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Spleen_weight.csv"
spleen <- read.csv(text = getURL(spleenURL))
names(spleen)[names(spleen)%in%"mouseID"] <- "EH_ID"
spleen <- merge(spleen, stab)
spleen$dpi <- as.numeric(gsub("dpi|dip", "", spleen$dpi.diss))

## weight --------------------------------------------------------
weightURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Weight_Expe_Jan_2017.csv"
weight <- read.csv(text = getURL(weightURL))

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
oocystsURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Clean_oocyst_data.csv"
oocysts <- read.csv(text = getURL(oocystsURL))

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
RtissueURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Eimeria-NMRI_Relative%20quantification_clean.csv"
Rtissue <- read.csv(text = getURL(RtissueURL))

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

histURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Histo.csv"
hist <- read.csv(text = getURL(histURL))
#add means and sums columns
hist$mMLS <- rowMeans(hist[, grepl("^MLS", colnames(hist))])
hist$sMLS <- rowSums(hist[, grepl("^MLS", colnames(hist))])

names(hist)[names(hist)%in%"Sample.ID"] <- "EH_ID"
all.data <- merge(all.data, hist, all=TRUE)

#order factors
all.data$dpi.diss <- factor(all.data$dpi.diss, levels = c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi"))

# test strain effect on IFNy CEWE and MES
IFNy.MES <- dplyr::select(all.data, EH_ID, dpi.diss, inf.strain, IFNy_MES)
IFNy.MES <- dplyr::distinct(IFNy.MES)
modIFNyMES <- lme4::lmer(IFNy_MES~inf.strain + (1|dpi.diss), data = IFNy.MES)
summary(modIFNyMES)

# graph
ggplot(data = all.data, aes(x = dpi.diss, y = IFNy_SPL, color = inf.strain, group = inf.strain)) + 
  geom_jitter(width = 0.2) +
  geom_smooth(se = F) +
  facet_grid(~inf.strain)

ggplot(data = all.data, aes(x = dpi.diss, y = IFNy_MES, color = inf.strain, group = dpi.diss & inf.strain)) + 
  geom_jitter(width=0.2) +
  geom_smooth(se = F) +
  facet_grid(~inf.strain)

# seems that spleen is not responsing to infection strain, MES better so let's explore that
# delta vs IFNy in MES 
ggplot(data = all.data, aes(x = PH.delta, y = IFNy_MES, color = inf.strain)) + 
  geom_jitter(width=0.2) +
  geom_smooth(se = F) +
  facet_grid(~inf.strain)
# oocyst counts vs IFNy in MES
ggplot(data = all.data, aes(x = Total.oocysts.g, y = IFNy_MES, color = inf.strain)) + 
  geom_jitter(width=0.2) +
  facet_grid(~inf.strain)
# that's some odd missing values

# weight vs IFNy MES
ggplot(data = all.data, aes(x = perc_of_dpi1, y = IFNy_MES, color = inf.strain)) + 
  geom_jitter(width=0.2) +
  facet_grid(~inf.strain)

# 
DT <- subset(all.data, all.data$Total.oocysts.g > 0 )

