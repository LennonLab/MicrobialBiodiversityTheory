rm(list = ls())
getwd()
setwd("~/github/MicroMETE/R/")

require(vegan)
library(reshape2)
library(plyr)
library(lme4)
library(ggplot2)
library(sjPlot)
library(VennDiagram)
# Import HMP metadata file
HMP_meta <- read.csv("../data/HMP-Data/ppAll_V35_map.txt",sep="\t")  # read csv file 
#HMP_meta[,c('NAP')]
# add suffix to NAP labels
HMP_meta$NAP <- sub("$", ".may1", HMP_meta$NAP )

# Import CSV of N, S, Nmax, r2, etc for HMP for each model tested
HMP_BS <- read.table("../data/NSR2/geom_HMP_NSR2.txt")
HMP_METE <- read.table("../data/NSR2/mete_HMP_NSR2.txt")
HMP_Zipf <- read.table("../data/NSR2/zipf_HMP_NSR2.txt")

####NOTE: RSID is the individual.
# add labels to columns
colnames(HMP_BS) <- c('site','N_BS','S_BS', 'Nmax_BS','R2_BS', 'NAP')
colnames(HMP_METE) <- c('site','N','S', 'Nmax','R2_METE', 'NAP')
colnames(HMP_Zipf) <- c('site','N','S', 'Nmax','gamma','R2_Zipf', 'NAP')
HMP_Zipf_keeps <- c('gamma','R2_Zipf', 'NAP')
HMP_Zipf_subset <- HMP_Zipf[HMP_Zipf_keeps]

# Select metadata columns we're interested in 
keeps <- c("NAP","Sex", "HMPBodySite", "HMPBodySubsite", "RSID", "VisitNo")
HMP_meta_subset <- HMP_meta[keeps]

# Merge into single matrix
data_merged <- merge(HMP_BS, HMP_METE, by = "NAP", all = TRUE)
data_merged <- merge(data_merged, HMP_Zipf_subset, by = "NAP", all = TRUE)
data_merged <- merge(HMP_meta_subset, data_merged)
drops <- c('N_BS','S_BS', 'Nmax_BS',"site.x" ,"site.y" )
data_merged <- data_merged[,!(names(data_merged) %in% drops)]
# Drop duplicate columns
data_merged_no_dups <- data_merged[!duplicated(data_merged[,1]),]

data_merged_no_dups_numeric <- data_merged_no_dups[,c("R2_BS","N","S", "Nmax", "R2_METE", 
                                                      "gamma", "R2_Zipf", "NAP", "RSID", "VisitNo",
                                                      "Sex", "HMPBodySite", "HMPBodySubsite")]


data_merged_no_dups_numeric$AvgAbund <- data_merged_no_dups_numeric$N/ data_merged_no_dups_numeric$S


log.cols <- c("N","AvgAbund", "S", "Nmax")
data_merged_no_dups_numeric[log.cols] <- log(data_merged_no_dups_numeric[log.cols], 10)
data <- data_merged_no_dups_numeric
data[, 'RSID'] <- as.factor(data[, 'RSID'])
data <- data[complete.cases(data),]


lm.geom <- lm(R2_BS ~ N * S, data = data)
lm.mete <- lm(R2_METE ~ N * S, data = data)
lm.zipf <- lm(R2_Zipf ~ N * S, data = data)


lmer.geom <- lmer(R2_BS ~ N * S +  (1|Sex/HMPBodySite/HMPBodySubsite), data = data, REML=FALSE)
lmer.mete <- lmer(R2_METE ~ N * S +  (1|Sex/HMPBodySite/HMPBodySubsite), data = data, REML=FALSE)
lmer.zipf <- lmer(R2_Zipf ~ N * S +  (1|Sex/HMPBodySite/HMPBodySubsite), data = data, REML=FALSE)



anova(lmer.geom, lm.geom)
summary(lmer.mete.noSex)
anova(lmer.mete, lm.mete)
summary(lmer.zipf.noSex)
anova(lmer.zipf, lm.zipf)


#lmer.geom.noSex <- lmer(R2_BS ~ N * S +  (1|HMPBodySite/HMPBodySubsite), data = data)
#lmer.mete.noSex <- lmer(R2_METE ~ N * S +  (1|HMPBodySite/HMPBodySubsite), data = data)
#lmer.zipf.noSex <- lmer(R2_Zipf ~ N * S +  (1|HMPBodySite/HMPBodySubsite), data = data)

#lmer.geom.noSex.dummy <- lmer(R2_BS ~ 1 + (1|HMPBodySite/HMPBodySubsite), data = data)
#lmer.mete.noSex.dummy <- lmer(R2_METE ~ 1 +  (1|HMPBodySite/HMPBodySubsite), data = data)
#lmer.zipf.noSex.dummy <- lmer(R2_Zipf ~ 1 +  (1|HMPBodySite/HMPBodySubsite), data = data)


lm1 = lm(R2_BS ~ N * S, data = data)
lm2 = lm(R2_BS ~ N + S, data = data)
anova(lm2, lm1)

summary(lmer.geom.noSex)
anova(lmer.geom.noSex, lm.geom)
summary(lmer.mete.noSex)
anova(lmer.mete.noSex, lm.mete)
summary(lmer.zipf.noSex)
anova(lmer.zipf.noSex, lm.zipf)


var.comp.mete <- ldply(VarCorr(lmer.mete.noSex))
var.residual.mete <-attr(VarCorr(lmer.mete.noSex), "sc")^2
names(var.comp.mete) <- c("Factor", "Variance")
var.comp.mete <- rbind(var.comp.mete, data.frame("Factor" = "Residual", "Variance" = attr(VarCorr(lmer.mete.noSex), "sc") ^ 2))
var.comp.mete$Percent <- round(var.comp.mete$Variance / sum(var.comp.mete$Variance) * 100, 1)
attr(var.comp.mete, "mer") <- lmer.mete.noSex


var.comp.zipf <- ldply(VarCorr(lmer.zipf.noSex))
var.residual.zipf <-attr(VarCorr(lmer.zipf.noSex), "sc")^2
names(var.comp.zipf) <- c("Factor", "Variance")
var.comp.zipf <- rbind(var.comp.zipf, data.frame("Factor" = "Residual", "Variance" = attr(VarCorr(lmer.zipf.noSex), "sc") ^ 2))
var.comp.zipf$Percent <- round(var.comp.zipf$Variance / sum(var.comp.zipf$Variance) * 100, 1)
attr(var.comp.zipf, "mer") <- lmer.zipf.noSex

var.comp.geom <- ldply(VarCorr(lmer.geom.noSex))
var.residual.geom <-attr(VarCorr(lmer.geom.noSex), "sc")^2
names(var.comp.geom) <- c("Factor", "Variance")
var.comp.geom <- rbind(var.comp.geom, data.frame("Factor" = "Residual", "Variance" = attr(VarCorr(lmer.geom.noSex), "sc") ^ 2))
var.comp.geom$Percent <- round(var.comp.geom$Variance / sum(var.comp.geom$Variance) * 100, 1)
attr(var.comp.geom, "mer") <- lmer.geom.noSex


grid.newpage()
draw.pairwise.venn(43.0, 1.6, 1.6, category = c("Body Site", "Body Subsite"), lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))


grid.newpage()
draw.pairwise.venn(36.0, 1.7, 1.7, category = c("Body Site", "Body Subsite"), lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))


grid.newpage()
draw.pairwise.venn(2.4, 2.4, 2.4, category = c("Body Site", "Body Subsite"), lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

