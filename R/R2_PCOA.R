rm(list = ls())
getwd()
setwd("~/github/MicroMETE/R/")

require("vegan")
library(reshape2)
library(plyr)


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
colnames(HMP_METE) <- c('site','N_METE','S_METE', 'Nmax_METE','R2_METE', 'NAP')
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




# Now do it for the Zipf dataset


# Now we can apply regression techniques.



# Now we can construct a covariance matrix
# merge two variables together 
#df.m<-melt(data_merged_no_nas_no_dups,id.vars=c("Sex","HMPBodySubsite"))

