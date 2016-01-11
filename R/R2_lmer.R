rm(list = ls())
getwd()
setwd("~/github/MicroMETE/R/")

require("vegan")
library(reshape2)
library(plyr)
library(lme4)

# Import HMP metadata file
HMP_meta <- read.csv("../data/HMP-Data/ppAll_V35_map_noTimeseries.txt",sep="\t")  # read csv file 
#HMP_meta[,c('NAP')]
# add suffix to NAP labels
HMP_meta$NAP <- sub("$", ".may1", HMP_meta$NAP )

# Import CSV of N, S, Nmax, r2, etc for HMP for each model tested
HMP_BS <- read.table("../data/NSR2/geom_HMP_NSR2.txt")
HMP_METE <- read.table("../data/NSR2/mete_HMP_NSR2.txt")
HMP_Zipf <- read.table("../data/NSR2/zipf_HMP_NSR2_subset.txt")

# add labels to columns
colnames(HMP_BS) <- c('site','N_BS','S_BS', 'Nmax_BS','R2_BS', 'NAP')
colnames(HMP_METE) <- c('site','N','S', 'Nmax','R2_METE', 'NAP')
colnames(HMP_Zipf) <- c('site','N','S', 'Nmax','gamma','R2_Zipf', 'NAP')
HMP_Zipf_keeps <- c('gamma','R2_Zipf', 'NAP')
HMP_Zipf_subset <- HMP_Zipf[HMP_Zipf_keeps]

# Select metadata columns we're interested in 
keeps <- c("NAP","Sex", "HMPBodySite", "HMPBodySubsite")
HMP_meta_subset <- HMP_meta[keeps]

# Merge into single matrix
data_merged <- merge(HMP_BS, HMP_METE, by = "NAP", all = TRUE)
data_merged <- merge(data_merged, HMP_Zipf_subset, by = "NAP", all = TRUE)
data_merged <- merge(HMP_meta_subset, data_merged)
drops <- c('N_BS','S_BS', 'Nmax_BS',"site.x" ,"site.y" )
data_merged <- data_merged[,!(names(data_merged) %in% drops)]
# Drop duplicate columns
data_merged_no_dups <- data_merged[!duplicated(data_merged[,1]),]

data_merged_no_dups_numeric <- data_merged_no_dups[,c("R2_BS","N","S", "Nmax", 
                                                      "R2_METE", "gamma", "R2_Zipf", "NAP",
                                                      "Sex", "HMPBodySite", "HMPBodySubsite")]


#  Log transform N0, S0, and Nmax
log.cols <- c("N","S", "Nmax")
data_merged_no_dups_numeric[log.cols] <- log(data_merged_no_dups_numeric[log.cols], 10)

data <- data_merged_no_dups_numeric

# now we can start. We're using the lme4 package, which allows for mixed effects.
# We're are considering sex, body site, and body sub-site as random effects
# We are considering N, S, and Nmax as fixed effects.

# The goal is to create three final models (one for each ecological model tested)

########
###Zipf#
########

# First, start by comparing N and S

lmN.zipf <- lm(R2_Zipf ~ N, data = data)
lmNS.zipf <- lm(R2_Zipf ~ N + S, data = data)
anova(lmN.zipf, lmNS.zipf)

# From the output we see that there's a significant difference, so keep S.
# Now test for an interaction (we honestly expect there to be one, but just to check...)
lmNS.int.zipf <- lm(R2_Zipf ~ N * S, data = data)
anova(lmNS.int.zipf, lmNS.zipf)

# Our expectations are correct. 
# Now add in Nmax
lmNSNmax.int.zipf <- lm(R2_Zipf ~ N * S + Nmax, data = data)
anova(lmNS.int.zipf, lmNSNmax.int.zipf)

# So keeep Nmax, now test for interaction. 
lmNSNmax.int.int.zipf <- lm(R2_Zipf ~ N * S * Nmax, data = data)
anova(lmNSNmax.int.zipf, lmNSNmax.int.int.zipf)

# Yeap, there's an interaction, but that's expected. 

# Now that that's compelted, we add in our random effects
# First, sex
lmer1.zipf <- lmer(R2_Zipf ~ N * S * Nmax +  (1|Sex), data = data)
summary(lmer1.zipf)
AIC(lmer1.zipf,lmNSNmax.int.zipf)
# A very small amount of the variance is explained at the Sex level
# AIC is also lower for initial model, so we can reject adding the random effect

# Let's look at Body site
lmer2.zipf <- lmer(R2_Zipf ~ N * S * Nmax +  (1|HMPBodySite), data = data)
summary(lmer2.zipf)
AIC(lmer2.zipf,lmNSNmax.int.int.zipf)

# So 4% of the variance is explained at the HMPBodySite level

# Now, examine sub-site
lmer3.zipf <- lmer(R2_Zipf ~ N * S * Nmax +  (1|HMPBodySubsite), data = data)
summary(lmer3.zipf)
AIC(lmer3.zipf,lmNSNmax.int.int.zipf)

# Body sub-site somewhat decreases the AIC, but this is also a multilevel dataset
# So we should nest sub-site within site
lmer4.zipf <- lmer(R2_Zipf ~ N * S * Nmax +  (1|HMPBodySite/HMPBodySubsite), data = data)
summary(lmer4.zipf)
AIC(lmer4.zipf,lmNSNmax.int.int.zipf)

# So AIC is slightly lower with Bodysite and sub-site included as random effects



########
###METE#
########


lmN.mete <- lm(R2_METE ~ N, data = data)
lmNS.mete <- lm(R2_METE ~ N + S, data = data)
anova(lmN.mete, lmNS.mete)

# From the output we see that there's a significant difference, so keep S.
# Now test for an interaction (we honestly expect there to be one, but just to check...)
lmNS.int.mete <- lm(R2_METE ~ N * S, data = data)
anova(lmNS.int.mete, lmNS.mete)

# Our expectations are correct. 
# Now add in Nmax
lmNSNmax.int.mete <- lm(R2_METE ~ N * S + Nmax, data = data)
anova(lmNS.int.mete, lmNSNmax.int.mete)

# So keeep Nmax, now test for interaction. 
lmNSNmax.int.int.mete <- lm(R2_METE ~ N * S * Nmax, data = data)
anova(lmNS.int.mete, lmNSNmax.int.int.mete)

# Yeap, there's an interaction, but that's expected. 

# Now that that's compelted, we add in our random effects
# First, sex
lmer1.mete <- lmer(R2_METE ~ N * S * Nmax +  (1|Sex), data = data)
summary(lmer1.mete)
AIC(lmer1.mete,lmNSNmax.int.int.mete)

# A very small amount of the variance is explained at the Sex level
# AIC is also lower for initial model, so we can reject adding the random effect

# Let's look at Body site
lmer2.mete <- lmer(R2_METE ~ N * S * Nmax +  (1|HMPBodySite), data = data)
summary(lmer2.mete)
AIC(lmer2.mete,lmNSNmax.int.int.mete)

# So ~9% of the variance is explained at the HMPBodySite level
# Also, the AIC is higher with that random effect, suggesting we should keep it

# Now, examine sub-site
lmer3.mete <- lmer(R2_METE ~ N * S * Nmax +  (1|HMPBodySubsite), data = data)
summary(lmer3.mete)
AIC(lmer3.mete,lmNSNmax.int.int.mete)

# So ~12% of the variance is explained at the HMPBodySite level
# Also, the AIC is higher with that random effect, suggesting we should keep it

# Body sub-site somewhat decreases the AIC, but this is also a multilevel dataset
# So we should nest sub-site within site
lmer4.mete <- lmer(R2_METE ~ N * S * Nmax +  (1|HMPBodySite/HMPBodySubsite), data = data)
summary(lmer4.mete)
AIC(lmer4.mete,lmNSNmax.int.int.mete)
# Subsite and body-site explain a lot of the vairance, so we should keep them. 


########
###BS###
########

lmN.BS <- lm(R2_BS ~ N, data = data)
lmNS.BS <- lm(R2_BS ~ N + S, data = data)
anova(lmN.BS, lmNS.BS)

# From the output we see that there's a significant difference, so keep S.
# Now test for an interaction (we honestly expect there to be one, but just to check...)
lmNS.int.BS <- lm(R2_BS ~ N * S, data = data)
anova(lmNS.int.BS, lmNS.BS)

# Our expectations are correct. 
# Now add in Nmax
lmNSNmax.int.BS <- lm(R2_BS ~ N * S + Nmax, data = data)
anova(lmNS.int.BS, lmNSNmax.int.BS)

# So keeep Nmax, now test for interaction. 
lmNSNmax.int.int.BS <- lm(R2_BS ~ N * S * Nmax, data = data)
anova(lmNS.int.BS, lmNSNmax.int.int.BS)

# Yeap, there's an interaction, but that's expected. 

# Now that that's compelted, we add in our random effects
# First, sex
lmer1.BS <- lmer(R2_BS ~ N * S * Nmax +  (1|Sex), data = data)
summary(lmer1.BS)
AIC(lmer1.BS,lmNSNmax.int.int.BS)
# A very small amount of the variance is explained at the Sex level
# AIC is also lower for initial model, so we can reject adding the random effect of Sex

# Let's look at Body site
lmer2.BS <- lmer(R2_BS ~ N * S * Nmax +  (1|HMPBodySite), data = data)
summary(lmer2.BS)
AIC(lmer2.BS,lmNSNmax.int.int.BS)

# So ~16% of the variance is explained at the HMPBodySite level
# Also, the AIC is higher with that random effect, suggesting we should keep it

# Now, examine sub-site
lmer3.BS <- lmer(R2_BS ~ N * S * Nmax +  (1|HMPBodySubsite), data = data)
summary(lmer3.BS)
AIC(lmer3.BS,lmNSNmax.int.int.BS)

# Body sub-site  decreases a lot the AIC, and ~27% variance is explained by sub-site
# but this is also a multilevel dataset, so we should nest sub-site within site
lmer4.BS <- lmer(R2_BS ~ N * S * Nmax +  (1|HMPBodySite/HMPBodySubsite), data = data)
summary(lmer4.BS)
AIC(lmer4.BS,lmNSNmax.int.int.BS)

# Subsite and body-site explain a lot of the vairance, so we should keep them. 

plot(lmer4.zipf )

# What to do next? 

# plot % variance explained? 

# AIC isn't the best criteria, 
# We could run a  parametric bootstrap using the posterior distribution of our model





