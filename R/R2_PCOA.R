rm(list = ls())
getwd()
setwd("~/github/MicroMETE/R/")

require(vegan)
library(reshape2)
library(plyr)
library(lme4)
library(ggplot2)
library(sjPlot)
# Import HMP metadata file
HMP_meta <- read.csv("../data/HMP-Data/ppAll_V35_map_noTimeseries.txt",sep="\t")  # read csv file 
#HMP_meta[,c('NAP')]
# add suffix to NAP labels
HMP_meta$NAP <- sub("$", ".may1", HMP_meta$NAP )

# Import CSV of N, S, Nmax, r2, etc for HMP for each model tested
HMP_BS <- read.table("../data/NSR2/geom_HMP_NSR2.txt")
HMP_METE <- read.table("../data/NSR2/mete_HMP_NSR2.txt")
HMP_Zipf <- read.table("../data/NSR2/zipf_HMP_NSR2.txt")

# add labels to columns
colnames(HMP_BS) <- c('site','N_BS','S_BS', 'Nmax_BS','R2_BS', 'NAP')
colnames(HMP_METE) <- c('site','N','S', 'Nmax','R2_METE', 'NAP')
colnames(HMP_Zipf) <- c('site','N','S', 'Nmax','gamma','R2_Zipf', 'NAP')
HMP_Zipf_keeps <- c('gamma','R2_Zipf', 'NAP')
HMP_Zipf_subset <- HMP_Zipf[HMP_Zipf_keeps]

# Select metadata columns we're interested in 
keeps <- c("NAP","Sex", "HMPBodySite", "HMPBodySubsite", "RSID")
HMP_meta_subset <- HMP_meta[keeps]

# Merge into single matrix
data_merged <- merge(HMP_BS, HMP_METE, by = "NAP", all = TRUE)
data_merged <- merge(data_merged, HMP_Zipf_subset, by = "NAP", all = TRUE)
data_merged <- merge(HMP_meta_subset, data_merged)
drops <- c('N_BS','S_BS', 'Nmax_BS',"site.x" ,"site.y" )
data_merged <- data_merged[,!(names(data_merged) %in% drops)]
# Drop duplicate columns
data_merged_no_dups <- data_merged[!duplicated(data_merged[,1]),]
data_merged_no_dups <-  data_merged_no_dups[!duplicated(data_merged_no_dups[,5]),]
data_merged_no_dups_numeric <- data_merged_no_dups[,c("R2_BS","N","S", "Nmax", 
                                                      "R2_METE", "gamma", "R2_Zipf", "NAP", "RSID",
                                                      "Sex", "HMPBodySite", "HMPBodySubsite")]


data_merged_no_dups_numeric$AvgAbund <- data_merged_no_dups_numeric$N/ data_merged_no_dups_numeric$S


log.cols <- c("N","AvgAbund", "S", "Nmax")
data_merged_no_dups_numeric[log.cols] <- log(data_merged_no_dups_numeric[log.cols], 10)
data <- data_merged_no_dups_numeric
data[, 'RSID'] <- as.factor(data[, 'RSID'])
data.nonas <- data[complete.cases(data),]
# start PCoA
#####################
####Broken-stick#####
#####################
r2.matrix <- data.nonas[,c("R2_BS","R2_METE", "R2_Zipf")]
r2.db <- vegdist(r2.matrix, method = "euclidean", upper = TRUE, diag = TRUE)
r2.pcoa <- cmdscale(r2.db, eig = TRUE, k = 3) 


explainvar1 <- round(r2.pcoa$eig[1] / sum(r2.pcoa$eig), 3) * 100
explainvar2 <- round(r2.pcoa$eig[2] / sum(r2.pcoa$eig), 3) * 100
explainvar3 <- round(r2.pcoa$eig[3] / sum(r2.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)

# Plot Eigenvalues
plot(r2.pcoa$eig, xlab = "PCoA Axis", ylab = "Eigenvalue", 
     las = 1, cex.lab = 1.5, pch = 16)

# Add Expectation based on Kaiser-Guttman criterion and Broken Stick Model
abline(h = mean(r2.pcoa$eig), lty = 2, lwd = 2, col = "blue")
b.stick <- bstick(29, sum(r2.pcoa$eig))
lines(1:29, b.stick, type = "l", lty = 4, lwd = 2, col = "red")

# Add Legend
legend("topright", legend = c("Avg Eigenvalue", "Broken-Stick"), 
       lty = c(2, 4), bty = "n", col = c("blue", "red"))

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)

# Initiate Plot
plot(r2.pcoa$points[ ,1], r2.pcoa$points[ ,2], ylim = c(-0.2, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
points(r2.pcoa$points[ ,1], r2.pcoa$points[ ,2],
       pch = 19, cex = 3, bg = "gray", col = "gray")
text(r2.pcoa$points[ ,1], r2.pcoa$points[ ,2], 
     labels = row.names(geom.pcoa$points))



# Define Environmental Matrix
r2.env <- data.nonas[,c("N","S", "Nmax", "AvgAbund", "RSID", "Sex", "HMPBodySite", "HMPBodySubsite")]
r2.env$Sex <- factor(r2.env$Sex)
r2.env$RSID <-  factor(r2.env$RSID)
r2.env$HMPBodySite <-  factor(r2.env$HMPBodySite)
r2.env$HMPBodySubsite <-  factor(r2.env$HMPBodySubsite)

r2.env <- data.matrix(r2.env)
# factor
# Conduct CCA 

doubs.rda <- rda (r2.matrix ~ r2.env)

#doubs.cca <- vegan::cca(r2.matrix ~ r2.env, scale = T)

# Permutation Tests
anova(rda.vasc, by = "axis")
rda.fit <- envfit(doubs.rda, r2.env, perm = 999)
#cca.fit

# Calculate Explained Variation
rda.explainvar1 <- round(doubs.rda$CCA$eig[1] / 
                           sum(c(doubs.rda$CCA$eig, doubs.rda$CA$eig)), 3) * 100
rda.explainvar2 <- round(doubs.rda$CCA$eig[2] / 
                           sum(c(doubs.rda$CCA$eig, doubs.rda$CA$eig)), 3) * 100

# Define Plot Parameters
par(mar = c(5, 5, 4, 4) + 0.1)

# Initiate Plot
plot(scores(doubs.rda, display = "wa"), xlim = c(-3.5, 2), ylim = c(-3.2, 3.2),
     xlab = paste("RDA 1 (", rda.explainvar1, "%)", sep = ""),
     ylab = paste("RDA 2 (", rda.explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
points(scores(doubs.rda, display = "wa"),
       pch = 19, cex = 3, bg = "gray", col = "gray")
#text(scores(doubs.rda, display = "wa"), 
#     labels = row.names(scores(doubs.rda, display = "wa")))

# Add Environmental Vectors
vectors <- scores(doubs.rda, display = "bp")
row.names(vectors) <- c("N","S", "Nmax", "RSID", "Sex", "HMPBodySite", "HMPBodySubsite")
arrows(0, 0, vectors[,1] * 2, vectors[, 2] * 2, 
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(vectors[,1] * 2, vectors[, 2] * 2, pos = 3, 
     labels = row.names(vectors))
axis(side = 3, lwd.ticks=2, cex.axis=1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[, 1])) * 2, labels = pretty(range(vectors[, 1])))
axis(side = 4, lwd.ticks=2, cex.axis=1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[, 2])) * 2, labels = pretty(range(vectors[, 2])))


